"""
Test spectra_wave_shape_2 for update_time
"""

import sys, os
import math, cmath
import numpy as np

import pytest

import shape_2
import corsys
import tfun

from spectral_wave_data import SpectralWaveData, SwdError, SwdFileCantOpenError, \
      SwdFileBinaryError, SwdFileDataError, SwdInputValueError, SwdAllocateError

assert sys.version_info > (3, 4)

# We check all permutations of...

dt_swd = 0.37        # Requested temporal spacing of SWD file (small due to polynomial variation)
nsteps = 20          # Requested number of time steps in SWD file.

impls = [1]          # List of spectral-wave-data implementations to test
dks = [0.3]          # delta(k) to test
ipols = [0, 1, 2]    # Temporal interpolation schemes
nswds = [3]          # number of non-zero spectral components in wave field
seeds = [0]          # Seed for random polynomial coefficients in temporal functions
#norders = [-1, 3]    # Taylor expansion order in z-direction
norders = [-1]    # Taylor expansion order in z-direction
ds = [3.0]           # Water depth
appsys = [(0.0, 0.0, 0.0, 0.0),  #'x0, y0, t0, beta'
          (0.0, 0.0, 1.23, 0.0)]

waves_param = []
for impl in impls:
    for dk in dks:
        for d in ds:
            for ipol in ipols:
                for nswd in nswds:
                    for seed in seeds:
                        for norder in norders:
                            for x0, y0, t0, beta in appsys:
                                waves_param.append({'impl':impl,
                                                    'dk':dk,
                                                    'd': d,
                                                    'ipol':ipol,
                                                    'nswd':nswd,
                                                    'seed':seed,
                                                    'norder':norder,
                                                    'x0':x0,
                                                    'y0':y0,
                                                    't0':t0,
                                                    'beta':beta})

def dicts2ids(ds):
    res = []
    for d in ds:
        s = ''
        for k, v in d.items():
            s += '_%s:%s' % (k,v)
        res.append(s)
    return res

@pytest.fixture(scope='function', params=waves_param, ids=dicts2ids(waves_param))
def make_waves(request, tmp_path_factory):
    # Setup code
    impl = request.param['impl']
    dk = request.param['dk']
    d = request.param['d']
    ipol = request.param['ipol']
    nswd = request.param['nswd']
    seed = request.param['seed']
    norder = request.param['norder']
    x0 = request.param['x0']
    y0 = request.param['y0']
    t0 = request.param['t0']
    beta = request.param['beta']
    mysys = corsys.CorSys(x0, y0, t0, beta)
    cfuns, hfuns = tfun.create_Tfuns_1d(nswd, ipol, seed)
    swd_anal = shape_2.Shape2(dk, nswd, d, cfuns, hfuns, norder, mysys)
    file_swd = 'impl_%i_dk_%.4f_nswd_%i_ipol_%i_norder_%i_x0_%.3f_y0_%.3f_t0_%.3f_beta_%.2f_seed_%i.swd' % \
               (impl, dk, nswd, ipol, norder, x0, y0, t0, beta, seed)
    swd_anal.tmpdir = str(tmp_path_factory.mktemp("swd"))
    file_swd = os.path.join(swd_anal.tmpdir, file_swd)
    file_swd_fail = os.path.join(swd_anal.tmpdir, file_swd[:-4] + '_FAIL.swd')
    swd_anal.write_swd(file_swd, dt=dt_swd, nsteps=nsteps)
    swd_anal.write_swd(file_swd_fail, dt=dt_swd, nsteps=nsteps, too_short_file=True)
    swd_num = SpectralWaveData(file_swd=file_swd, x0=x0, y0=y0, t0=t0, beta=beta,
                               rho=1025.0, nsumx=-1, nsumy=-1, impl=impl, ipol=ipol,
                               norder=norder, dc_bias=True)
    swd_num_fail = SpectralWaveData(file_swd=file_swd_fail, x0=x0, y0=y0, t0=t0, beta=beta,
                                    rho=1025.0, nsumx=-1, nsumy=-1, impl=impl, ipol=ipol,
                                    norder=norder, dc_bias=True)
    # Inject class variables
    yield swd_anal, swd_num, swd_num_fail
    # Teardown code
    swd_num.close()
    swd_num_fail.close()


def test_waves(make_waves):
    swd_anal, swd_num, swd_num_fail = make_waves

    x = 2.0  # application x-positions for evaluation
    y = 3.0  # application y-positions for evaluation
    z = -4.0  # application z-positions for evaluation

    t0 = swd_num['t0']
    tmax = swd_num['tmax']
    assert math.isclose(tmax, dt_swd * (nsteps - 1) - t0, rel_tol=1e-04, abs_tol=1e-04)

    def checks(ts, analytic=True):
        res_phi = []
        res_elev = []
        for t in ts:
            swd_num.update_time(t)
            phi_num = swd_num.phi(x, y, z)
            if analytic:
                phi_anal = swd_anal.phi(x, y, z, t)
                assert math.isclose(phi_num, phi_anal, rel_tol=1e-03, abs_tol=1e-03)
            res_phi.append(phi_num)

            elev_num = swd_num.elev(x, y)
            if analytic:
                elev_anal = swd_anal.elev(x, y, t)
                assert math.isclose(elev_num, elev_anal, rel_tol=1e-03, abs_tol=1e-03)
            res_elev.append(elev_num)
        return np.array(res_phi), np.array(res_elev)

    # Test time series with different sampling density compared to dt_swd.
    # Then we repeat the calculations with random shuffling of time
    # instances to check if we get the very same values. The random shuffle will test
    # effects of positive and negative time stepping of small and large magnitudes.

    dt_div_dt_swd = [0.34, 0.44, 0.22, 0.6, 0.0, 1.5, 0.3, 2.3, 3.1, 4.0, 5.0, 0.2]
    ts = [0.0]
    for f in dt_div_dt_swd:
        ts.append(ts[-1] + f * dt_swd)
    ts += [tmax - 2.3 * dt_swd, tmax - 0.3 * dt_swd, tmax]
    if t0 > 0.0:
        ts = [t for t in ts if t <= tmax]
    ts = np.array(ts)
    nsamp = ts.size
    if swd_num['ipol'] < 2:
        res_phi, res_elev = checks(ts, analytic=True)
    else:
        # ipol=2 is not really consistent with a global polynomial.
        res_phi, res_elev = checks(ts, analytic=False)
    for j in range(2):
        ip = np.random.permutation(nsamp)
        ts_2 = ts[ip]
        res_phi_2, res_elev_2 = checks(ts_2, analytic=False)
        for i in ip:
            assert math.isclose(res_phi_2[i], res_phi[ip[i]], rel_tol=1e-10, abs_tol=1e-10)
            assert math.isclose(res_elev_2[i], res_elev[ip[i]], rel_tol=1e-10, abs_tol=1e-10)

    with pytest.raises(SwdInputValueError):
        t = - 0.01 * dt_swd - t0  # Negative t_swd
        swd_num.update_time(t)

    with pytest.raises(SwdInputValueError):
        t = tmax + 0.01*dt_swd
        swd_num.update_time(t)

    with pytest.raises(SwdFileDataError):
        t = tmax - dt_swd   # SWD file is too short. (Ended prematurely)
        swd_num_fail.update_time(t)
