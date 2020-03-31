"""
Test spectra_wave_shape_2 against analytical fields evaluated by sympy
"""

import sys, os
import math, cmath

import pytest

import shape_2
import corsys
import tfun

from spectral_wave_data import SpectralWaveData

assert sys.version_info > (3, 4)

# We check all permutations of...

impls = [1]          # List of spectral-wave-data implementations to test
dks = [0.3]          # delta(k) to test
ipols = [0, 1, 2]    # Temporal interpolation schemes
nswds = [3]          # number of non-zero spectral components in wave field
seeds = [0]          # Seed for random polynomial coefficients in temporal functions
#norders = [-1, 3]   # Taylor expansion order in z-direction
norders = [-1]       # Taylor expansion order in z-direction
ds = [3.0]           # Water depth
appsys = [(0.0, 0.0, 0.0, 0.0),  #'x0, y0, t0, beta'
          (13.8, 0.0, 0.0, 0.0),
          (0.0, 18.2, 0.0, 0.0),
          (0.0, 0.0, 0.34, 0.0),
          (-5.3, 18.2, 0.0, 37.4),
          (-5.3, 18.2, 0.0, 137.4),
          (-5.3, 18.2, 0.0, 237.4),
          (-5.3, 18.2, 0.0, 337.4)]  # Coordinate transformations

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
                                                    'd':d,
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
    file_swd = 'impl_%i_dk_%.4f_d_%.4f_nswd_%i_ipol_%i_norder_%i_x0_%.3f_y0_%.3f_t0_%.3f_beta_%.2f_seed_%i.swd' % \
               (impl, dk, d, nswd, ipol, norder, x0, y0, t0, beta, seed)
    swd_anal.tmpdir = str(tmp_path_factory.mktemp("swd"))
    file_swd = os.path.join(swd_anal.tmpdir, file_swd)
    #print('file_swd = ', file_swd)
    #swd_anal.dump_spectral_fun(j=0, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=1, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=2, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=3, dt=0.01, tmax=1.0)
    swd_anal.write_swd(file_swd, dt=0.1, nsteps=11)
    swd_num = SpectralWaveData(file_swd, x0, y0, t0, beta, rho=1025.0, impl=impl, dc_bias=True)
    # Inject class variables
    #request.cls.swd_num = swd_num
    #request.cls.swd_anal = swd_anal
    yield swd_anal, swd_num
    # Teardown code
    #print('CLOSE: ', file_swd)
    swd_num.close()

xs = [0.0, 8.4]          # application x-positions for evaluation
ys = [0.0, -5.2]         # application y-positions for evaluation
zs = [0.0, -2.3, 1.7]    # application z-positions for evaluation
ts = [0.2, 0.23, 0.52]   # application time: First exact on 2dt, Then < dt, then > dt

def test_waves(make_waves):
    swd_anal, swd_num = make_waves

    for t in ts:

        swd_num.update_time(t)

        for x in xs:
            for y in ys:
                for z in zs:

                    phi_anal = swd_anal.phi(x, y, z, t)
                    phi_num = swd_num.phi(x, y, z)
                    assert math.isclose(phi_num, phi_anal, rel_tol=1e-04, abs_tol=1e-04)

                    phi_t_anal = swd_anal.phi_t(x, y, z, t)
                    phi_t_num = swd_num.phi_t(x, y, z)
                    assert math.isclose(phi_t_num, phi_t_anal, rel_tol=1e-04, abs_tol=1e-04)

                    stream_anal = swd_anal.stream(x, y, z, t)
                    stream_num = swd_num.stream(x, y, z)
                    assert math.isclose(stream_num, stream_anal, rel_tol=1e-04, abs_tol=1e-04)

                    grad_phi_anal = swd_anal.grad_phi(x, y, z, t)
                    grad_phi_num = swd_num.grad_phi(x, y, z)
                    assert math.isclose(grad_phi_num.x, grad_phi_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num.y, grad_phi_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num.z, grad_phi_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_anal = swd_anal.grad_phi_2nd(x, y, z, t)
                    g2nd_num = swd_num.grad_phi_2nd(x, y, z)
                    assert math.isclose(g2nd_num.xx, g2nd_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num.xy, g2nd_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num.xz, g2nd_anal['xz'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num.yy, g2nd_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num.yz, g2nd_anal['yz'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num.zz, g2nd_anal['zz'], rel_tol=1e-04, abs_tol=1e-04)

                    acc_euler_anal = swd_anal.acc_euler(x, y, z, t)
                    acc_euler_num = swd_num.acc_euler(x, y, z)
                    assert math.isclose(acc_euler_num.x, acc_euler_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_euler_num.y, acc_euler_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_euler_num.z, acc_euler_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    acc_particle_anal = swd_anal.acc_particle(x, y, z, t)
                    acc_particle_num = swd_num.acc_particle(x, y, z)
                    assert math.isclose(acc_particle_num.x, acc_particle_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_particle_num.y, acc_particle_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_particle_num.z, acc_particle_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    prs_anal = swd_anal.pressure(x, y, z, t, rho=1025.0, grav=9.81)
                    prs_num = swd_num.pressure(x, y, z)
                    assert math.isclose(prs_num, prs_anal, rel_tol=1e-03)

                    elev_anal = swd_anal.elev(x, y, t)
                    elev_num = swd_num.elev(x, y)
                    assert math.isclose(elev_num, elev_anal, rel_tol=1e-04, abs_tol=1e-04)

                    elev_t_anal = swd_anal.elev_t(x, y, t)
                    elev_t_num = swd_num.elev_t(x, y)
                    assert math.isclose(elev_t_num, elev_t_anal, rel_tol=1e-04, abs_tol=1e-04)

                    grad_elev_anal = swd_anal.grad_elev(x, y, t)
                    grad_elev_num = swd_num.grad_elev(x, y)
                    assert math.isclose(grad_elev_num.x, grad_elev_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_elev_num.y, grad_elev_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_elev_num.z, grad_elev_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_elev_anal = swd_anal.grad_elev_2nd(x, y, t)
                    g2nd_elev_num = swd_num.grad_elev_2nd(x, y)
                    assert math.isclose(g2nd_elev_num.xx, g2nd_elev_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_elev_num.xy, g2nd_elev_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_elev_num.yy, g2nd_elev_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)

                    if abs(swd_anal.sys.t0 - 0.34) < 1.0e-6 and abs(t - 0.52) < 1.0e-6:
                        # Check if strip method produce a swd file with consistent kinematic results.
                        file_swd_strip = os.path.join(swd_anal.tmpdir, 'stripped.swd')
                        swd_num.strip(tmin=0.31, tmax=0.61, file_swd=file_swd_strip)
                        swd_num_strip = SpectralWaveData(file_swd_strip,
                                                         x0=swd_anal.sys.x0,
                                                         y0=swd_anal.sys.y0,
                                                         t0=0.0,
                                                         beta=swd_anal.sys.beta,
                                                         rho=1025.0, impl=impl, dc_bias=True)
                        nstrip = swd_num_strip.get('nstrip')
                        assert nstrip > 0

                        t_strip = (swd_anal.sys.t0 + t) - nstrip * 0.1 # dt=0.1 in swd file, nstrip includes t=0.0
                        swd_num_strip.update_time(t_strip)

                        phi_num_strip = swd_num_strip.phi(x, y, z)
                        assert math.isclose(phi_num, phi_num_strip, rel_tol=1e-04, abs_tol=1e-04)

                        elev_t_num_strip = swd_num_strip.elev_t(x, y)
                        assert math.isclose(elev_t_num, elev_t_num_strip, rel_tol=1e-04, abs_tol=1e-04)

                        swd_num_strip.close()

                    file_convergence = os.path.join(swd_anal.tmpdir, 'file_convergence.csv')
                    swd_num.convergence(x, y, z, csv=file_convergence)
                    with open(file_convergence, 'r') as f:
                        lines = f.read().splitlines()
                        words = lines[-1].split(',')

                    phi_x_convergence = float(words[1])
                    phi_y_convergence = float(words[2])
                    phi_z_convergence = float(words[3])
                    assert math.isclose(grad_phi_num.x, phi_x_convergence, rel_tol=1e-07, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num.y, phi_y_convergence, rel_tol=1e-07, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num.z, phi_z_convergence, rel_tol=1e-07, abs_tol=1e-04)

                    elev_convergence = float(words[4])
                    assert math.isclose(elev_num, elev_convergence, rel_tol=1e-07, abs_tol=1e-04)

                    prs_convergence = float(words[5])
                    assert math.isclose(prs_num, prs_convergence, rel_tol=1e-07, abs_tol=1e-04)

        z_floor = swd_num.bathymetry(x=3.4, y=-5.3)
        assert math.isclose(z_floor, swd_anal.d, rel_tol=1e-5)

        nvec_floor = swd_num.bathymetry_nvec(x=3.4, y=-5.3)
        assert math.isclose(nvec_floor.x, 0.0, abs_tol=1e-5)
        assert math.isclose(nvec_floor.y, 0.0, abs_tol=1e-5)
        assert math.isclose(nvec_floor.z, 1.0, rel_tol=1e-5)

        nswd_file = swd_num.get('n')
        assert isinstance(nswd_file, int)
        assert nswd_file == 3

        dt_swd = swd_num.get('dt')
        assert isinstance(dt_swd, float)
        assert math.isclose(dt_swd, 0.1, rel_tol=1e-5)

        cid = swd_num.get('cid')
        assert isinstance(cid, str)
        assert cid == "{'ole':1, 'dole':2, 'doffen':3}"

