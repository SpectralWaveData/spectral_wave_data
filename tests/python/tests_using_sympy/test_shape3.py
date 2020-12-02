"""
Test spectra_wave_shape_2 against analytical fields evaluated by sympy
"""

import sys, os
import math, cmath

import pytest
import numpy as np

import shape_1
import shape_2
import shape_3
import corsys
import tfun

from spectral_wave_data import SpectralWaveData

assert sys.version_info > (3, 4)

# We check all permutations of...

impls = [1]          # List of spectral-wave-data implementations to test
dks = [0.3]          # delta(k) to test
ipols = [0, 1, 2]    # Temporal interpolation schemes
nswds = [3]          # number of non-zero spectral components in wave field
nsfs = [0, 1, 3]     # Number of sea floor offset points
seeds = [0]          # Seed for random polynomial coefficients in temporal functions
#norders = [-1, 3]   # Taylor expansion order in z-direction
norders = [-1]       # Taylor expansion order in z-direction
appsys = [(0.0, 0.0, 0.0, 0.0),  #'x0, y0, t0, beta'
          (13.8, 0.0, 0.0, 0.0),
          (0.0, 18.2, 0.0, 0.0),
          (0.0, 0.0, 0.34, 0.0),
          (-5.3, 18.2, 0.0, 37.4)]  # Coordinate transformations

waves_param = []
for impl in impls:
    for dk in dks:
        for ipol in ipols:
            for nsf in nsfs:
                for nswd in nswds:
                    for seed in seeds:
                        for norder in norders:
                            for x0, y0, t0, beta in appsys:
                                waves_param.append({'impl':impl,
                                                    'dk':dk,
                                                    'nsf':nsf,
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
    nsf = request.param['nsf']
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
    isf = 0
    if nsf == 0:
        nh = -1
        chfuns = []
        xsf = []
        zsf = []
        swd_anal = shape_3.Shape3(dk, nswd, nh, isf, nsf, xsf, zsf, cfuns, chfuns, hfuns, norder, mysys)
    elif nsf == 1:
        nh = nswd
        chfuns = [] # Will be constructed based on flat bottom requirement
        xsf = [0.0]
        zsf = [-3.0]  # Water depth = 3m
        swd_anal = shape_3.Shape3(dk, nswd, nh, isf, nsf, xsf, zsf, cfuns, chfuns, hfuns, norder, mysys)
    else:
        nh = nswd - 1 # To check difference in loop counters
        chfuns, __ = tfun.create_Tfuns_1d(nh, ipol, seed + 1)
        xsf = np.linspace(0.0, 2*np.pi/dk, nsf)
        zsf = [-4.0 - 0.8 * i**2 for i in range(nsf)]
        swd_anal = shape_3.Shape3(dk, nswd, nh, isf, nsf, xsf, zsf, cfuns, chfuns, hfuns, norder, mysys)
    swd_anal.tmpdir = str(tmp_path_factory.mktemp("swd"))
    file_swd1 = os.path.join(swd_anal.tmpdir, "shp1.swd")
    file_swd2 = os.path.join(swd_anal.tmpdir, "shp2.swd")
    file_swd3 = os.path.join(swd_anal.tmpdir, "shp3.swd")

    if nsf == 0:
        swd_anal.write_swd1(file_swd1, dt=0.1, nsteps=11)
        swd_anal.check_swd_meta(file_swd1, 1, nswd, nh)
        swd_num1 = SpectralWaveData(file_swd1, x0, y0, t0, beta, rho=1025.0, impl=impl, dc_bias=True)
    else:
        swd_num1 = None
    if nsf == 1:
        swd_anal.write_swd2(file_swd2, dt=0.1, nsteps=11)
        swd_anal.check_swd_meta(file_swd2, 2, nswd, nh)
        swd_num2 = SpectralWaveData(file_swd2, x0, y0, t0, beta, rho=1025.0, impl=impl, dc_bias=True)
    else:
        swd_num2 = None

    swd_anal.write_swd3(file_swd3, dt=0.1, nsteps=11)
    swd_anal.check_swd_meta(file_swd3, 3, nswd, nh)
    swd_num3 = SpectralWaveData(file_swd3, x0, y0, t0, beta, rho=1025.0, impl=impl, dc_bias=True)

    yield swd_anal, swd_num1, swd_num2, swd_num3, nsf

    # Teardown code
    if nsf == 0:
        #print("CLOSE: ", file_swd1)
        swd_num1.close()
    if nsf == 1:
        #print("CLOSE: ", file_swd2)
        swd_num2.close()
    #print("CLOSE: ", file_swd3)
    swd_num3.close()

xs = [0.0, 8.4]          # application x-positions for evaluation
ys = [0.0, -5.2]         # application y-positions for evaluation
zs = [0.0, -2.3, 1.7]    # application z-positions for evaluation
ts = [0.2, 0.23, 0.52]   # application time: First exact on 2dt, Then < dt, then > dt

def test_waves(make_waves):
    swd_anal, swd_num1, swd_num2, swd_num3, nsf = make_waves

    #swd_anal.dump_spectral_fun(j=0, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=1, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=2, dt=0.01, tmax=1.0)
    #swd_anal.dump_spectral_fun(j=3, dt=0.01, tmax=1.0)

    for t in ts:

        if swd_num1 is not None:
            swd_num1.update_time(t)
        if swd_num2 is not None:
            swd_num2.update_time(t)
        swd_num3.update_time(t)

        for x in xs:
            for y in ys:
                for z in zs:

                    phi_anal = swd_anal.phi(x, y, z, t)
                    phi_num3 = swd_num3.phi(x, y, z)
                    assert math.isclose(phi_num3, phi_anal, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        phi_num1 = swd_num1.phi(x, y, z)
                        assert math.isclose(phi_num1, phi_num3, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        phi_num2 = swd_num2.phi(x, y, z)
                        assert math.isclose(phi_num2, phi_num3, rel_tol=1e-04, abs_tol=1e-04)

                    stream_anal = swd_anal.stream(x, y, z, t)
                    stream_num3 = swd_num3.stream(x, y, z)
                    assert math.isclose(stream_num3, stream_anal, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        stream_num1 = swd_num1.stream(x, y, z)
                        assert math.isclose(stream_num1, stream_num3, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        stream_num2 = swd_num2.stream(x, y, z)
                        assert math.isclose(stream_num2, stream_num3, rel_tol=1e-04, abs_tol=1e-04)

                    phi_t_anal = swd_anal.phi_t(x, y, z, t)
                    phi_t_num3 = swd_num3.phi_t(x, y, z)
                    assert math.isclose(phi_t_num3, phi_t_anal, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        phi_t_num1 = swd_num1.phi_t(x, y, z)
                        assert math.isclose(phi_t_num1, phi_t_num3, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        phi_t_num2 = swd_num2.phi_t(x, y, z)
                        assert math.isclose(phi_t_num2, phi_t_num3, rel_tol=1e-04, abs_tol=1e-04)

                    grad_phi_anal = swd_anal.grad_phi(x, y, z, t)
                    grad_phi_num3 = swd_num3.grad_phi(x, y, z)
                    assert math.isclose(grad_phi_num3.x, grad_phi_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num3.y, grad_phi_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num3.z, grad_phi_anal['z'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        grad_phi_num1 = swd_num1.grad_phi(x, y, z)
                        assert math.isclose(grad_phi_num1.x, grad_phi_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num1.y, grad_phi_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num1.z, grad_phi_num3.z, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        grad_phi_num2 = swd_num2.grad_phi(x, y, z)
                        assert math.isclose(grad_phi_num2.x, grad_phi_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num2.y, grad_phi_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num2.z, grad_phi_num3.z, rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_anal = swd_anal.grad_phi_2nd(x, y, z, t)
                    g2nd_num3 = swd_num3.grad_phi_2nd(x, y, z)
                    assert math.isclose(g2nd_num3.xx, g2nd_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num3.xy, g2nd_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num3.xz, g2nd_anal['xz'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num3.yy, g2nd_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num3.yz, g2nd_anal['yz'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_num3.zz, g2nd_anal['zz'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        g2nd_num1 = swd_num1.grad_phi_2nd(x, y, z)
                        assert math.isclose(g2nd_num1.xx, g2nd_num3.xx, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num1.xy, g2nd_num3.xy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num1.xz, g2nd_num3.xz, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num1.yy, g2nd_num3.yy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num1.yz, g2nd_num3.yz, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num1.zz, g2nd_num3.zz, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        g2nd_num2 = swd_num2.grad_phi_2nd(x, y, z)
                        assert math.isclose(g2nd_num2.xx, g2nd_num3.xx, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num2.xy, g2nd_num3.xy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num2.xz, g2nd_num3.xz, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num2.yy, g2nd_num3.yy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num2.yz, g2nd_num3.yz, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num2.zz, g2nd_num3.zz, rel_tol=1e-04, abs_tol=1e-04)

                    acc_euler_anal = swd_anal.acc_euler(x, y, z, t)
                    acc_euler_num3 = swd_num3.acc_euler(x, y, z)
                    assert math.isclose(acc_euler_num3.x, acc_euler_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_euler_num3.y, acc_euler_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_euler_num3.z, acc_euler_anal['z'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        acc_euler_num1 = swd_num1.acc_euler(x, y, z)
                        assert math.isclose(acc_euler_num1.x, acc_euler_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_euler_num1.y, acc_euler_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_euler_num1.z, acc_euler_num3.z, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        acc_euler_num2 = swd_num2.acc_euler(x, y, z)
                        assert math.isclose(acc_euler_num2.x, acc_euler_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_euler_num2.y, acc_euler_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_euler_num2.z, acc_euler_num3.z, rel_tol=1e-04, abs_tol=1e-04)

                    acc_particle_anal = swd_anal.acc_particle(x, y, z, t)
                    acc_particle_num3 = swd_num3.acc_particle(x, y, z)
                    assert math.isclose(acc_particle_num3.x, acc_particle_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_particle_num3.y, acc_particle_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(acc_particle_num3.z, acc_particle_anal['z'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        acc_particle_num1 = swd_num1.acc_particle(x, y, z)
                        assert math.isclose(acc_particle_num1.x, acc_particle_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_particle_num1.y, acc_particle_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_particle_num1.z, acc_particle_num3.z, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        acc_particle_num2 = swd_num2.acc_particle(x, y, z)
                        assert math.isclose(acc_particle_num2.x, acc_particle_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_particle_num2.y, acc_particle_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(acc_particle_num2.z, acc_particle_num3.z, rel_tol=1e-04, abs_tol=1e-04)

                    prs_anal = swd_anal.pressure(x, y, z, t, rho=1025.0, grav=9.81)
                    prs_num3 = swd_num3.pressure(x, y, z)
                    assert math.isclose(prs_num3, prs_anal, rel_tol=1e-03)
                    if swd_num1 is not None:
                        prs_num1 = swd_num1.pressure(x, y, z)
                        assert math.isclose(prs_num1, prs_num3, rel_tol=1e-03)
                    if swd_num2 is not None:
                        prs_num2 = swd_num2.pressure(x, y, z)
                        assert math.isclose(prs_num2, prs_num3, rel_tol=1e-03)

                    elev_anal = swd_anal.elev(x, y, t)
                    elev_num3 = swd_num3.elev(x, y)
                    assert math.isclose(elev_num3, elev_anal, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        elev_num1 = swd_num1.elev(x, y)
                        assert math.isclose(elev_num1, elev_num3, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        elev_num2 = swd_num2.elev(x, y)
                        assert math.isclose(elev_num2, elev_num3, rel_tol=1e-04, abs_tol=1e-04)

                    elev_t_anal = swd_anal.elev_t(x, y, t)
                    elev_t_num3 = swd_num3.elev_t(x, y)
                    assert math.isclose(elev_t_num3, elev_t_anal, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        elev_t_num1 = swd_num1.elev_t(x, y)
                        assert math.isclose(elev_t_num1, elev_t_num3, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        elev_t_num2 = swd_num2.elev_t(x, y)
                        assert math.isclose(elev_t_num2, elev_t_num3, rel_tol=1e-04, abs_tol=1e-04)

                    grad_elev_anal = swd_anal.grad_elev(x, y, t)
                    grad_elev_num3 = swd_num3.grad_elev(x, y)
                    assert math.isclose(grad_elev_num3.x, grad_elev_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_elev_num3.y, grad_elev_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_elev_num3.z, grad_elev_anal['z'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        grad_elev_num1 = swd_num1.grad_elev(x, y)
                        assert math.isclose(grad_elev_num1.x, grad_elev_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num1.y, grad_elev_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num1.z, grad_elev_num3.z, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        grad_elev_num2 = swd_num2.grad_elev(x, y)
                        assert math.isclose(grad_elev_num2.x, grad_elev_num3.x, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num2.y, grad_elev_num3.y, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num2.z, grad_elev_num3.z, rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_elev_anal = swd_anal.grad_elev_2nd(x, y, t)
                    g2nd_elev_num3 = swd_num3.grad_elev_2nd(x, y)
                    assert math.isclose(g2nd_elev_num3.xx, g2nd_elev_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_elev_num3.xy, g2nd_elev_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(g2nd_elev_num3.yy, g2nd_elev_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num1 is not None:
                        g2nd_elev_num1 = swd_num1.grad_elev_2nd(x, y)
                        assert math.isclose(g2nd_elev_num1.xx, g2nd_elev_num3.xx, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num1.xy, g2nd_elev_num3.xy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num1.yy, g2nd_elev_num3.yy, rel_tol=1e-04, abs_tol=1e-04)
                    if swd_num2 is not None:
                        g2nd_elev_num2 = swd_num2.grad_elev_2nd(x, y)
                        assert math.isclose(g2nd_elev_num2.xx, g2nd_elev_num3.xx, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num2.xy, g2nd_elev_num3.xy, rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num2.yy, g2nd_elev_num3.yy, rel_tol=1e-04, abs_tol=1e-04)

                    if abs(swd_anal.sys.t0 - 0.34) < 1.0e-6 and abs(t - 0.52) < 1.0e-6:
                        # Check if strip method produce a swd file with consistent kinematic results.
                        file_swd_strip = os.path.join(swd_anal.tmpdir, 'stripped.swd')
                        swd_num3.strip(tmin=0.31, tmax=0.61, file_swd=file_swd_strip)
                        swd_num3_strip = SpectralWaveData(file_swd_strip,
                                                          x0=swd_anal.sys.x0,
                                                          y0=swd_anal.sys.y0,
                                                          t0=0.0,
                                                          beta=swd_anal.sys.beta,
                                                          rho=1025.0, impl=impl, dc_bias=True)
                        nstrip = swd_num3_strip.get('nstrip')
                        assert nstrip > 0

                        t_strip = (swd_anal.sys.t0 + t) - nstrip * 0.1 # dt=0.1 in swd file, nstrip includes t=0.0
                        swd_num3_strip.update_time(t_strip)

                        phi_num3_strip = swd_num3_strip.phi(x, y, z)
                        assert math.isclose(phi_num3, phi_num3_strip, rel_tol=1e-04, abs_tol=1e-04)

                        elev_t_num3_strip = swd_num3_strip.elev_t(x, y)
                        assert math.isclose(elev_t_num3, elev_t_num3_strip, rel_tol=1e-04, abs_tol=1e-04)

                        swd_num3_strip.close()

                    file_convergence = os.path.join(swd_anal.tmpdir, 'file_convergence.csv')
                    swd_num3.convergence(x, y, z, csv=file_convergence)
                    with open(file_convergence, 'r') as f:
                        lines = f.read().splitlines()
                        words = lines[-1].split(',')

                    phi_x_convergence = float(words[1])
                    phi_y_convergence = float(words[2])
                    phi_z_convergence = float(words[3])
                    assert math.isclose(grad_phi_num3.x, phi_x_convergence, rel_tol=1e-07, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num3.y, phi_y_convergence, rel_tol=1e-07, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num3.z, phi_z_convergence, rel_tol=1e-07, abs_tol=1e-04)

                    elev_convergence = float(words[4])
                    assert math.isclose(elev_num3, elev_convergence, rel_tol=1e-07, abs_tol=1e-04)

                    prs_convergence = float(words[5])
                    assert math.isclose(prs_num3, prs_convergence, rel_tol=1e-07, abs_tol=1e-04)

        depth_anal = swd_anal.bathymetry(x_app=3.4, y_app=-5.3)
        depth_num3 = swd_num3.bathymetry(x=3.4, y=-5.3)
        if depth_anal < 0:
            assert depth_num3 < 0
        else:
            assert math.isclose(depth_num3, depth_anal, rel_tol=1e-5)

        nvec_floor_anal = swd_anal.bathymetry_nvec(x_app=3.4, y_app=-5.3)
        nvec_floor_num3 = swd_num3.bathymetry_nvec(x=3.4, y=-5.3)
        assert math.isclose(nvec_floor_anal['x'], nvec_floor_num3.x, abs_tol=1e-5)
        assert math.isclose(nvec_floor_anal['y'], nvec_floor_num3.y, abs_tol=1e-5)
        assert math.isclose(nvec_floor_anal['z'], nvec_floor_num3.z, rel_tol=1e-5)

        nswd_file = swd_num3.get('n')
        assert isinstance(nswd_file, int)
        assert nswd_file == 3

        dt_swd = swd_num3.get('dt')
        assert isinstance(dt_swd, float)
        assert math.isclose(dt_swd, 0.1, rel_tol=1e-5)

        dc_bias = swd_num3["dc_bias"]
        assert isinstance(dc_bias, bool)
        assert dc_bias is True

        cid = swd_num3.get('cid')
        assert isinstance(cid, str)
        assert cid == "{'ole':1, 'dole':2, 'doffen':3}"

