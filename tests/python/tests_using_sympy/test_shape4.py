"""
Test spectra_wave_shape_4 against analytical fields evaluated by sympy
"""

import sys, os, shutil
import math, cmath

import pytest
import numpy as np

import shape_4
import corsys
import tfun

from spectral_wave_data import SpectralWaveData

assert sys.version_info > (3, 4)

# We check all permutations of...

ipols = [0, 1, 2]    # Temporal interpolation schemes
seeds = [0]          # Seed for random polynomial coefficients in temporal functions
#norders = [-1, 3]   # Taylor expansion order in z-direction
norders = [-1]       # Taylor expansion order in z-direction
resolution = [#'nx, ny, dkx, dky'
              (3, 3, 0.3, 0.3),  # Symmetric space (applicable for impl=2 too)
              (3, 2, 0.3, 0.2)]  # General non-symmetric.
appsys = [(0.0, 0.0, 0.0, 0.0),  #'x0, y0, t0, beta'
          (13.8, 0.0, 0.0, 0.0),
          (0.0, 18.2, 0.0, 0.0),
          (0.0, 0.0, 0.34, 0.0),
          (-5.3, 18.2, 0.0, 37.4)]  # Coordinate transformations

waves_param = []
for nx, ny, dkx, dky in resolution:
    for ipol in ipols:
            for seed in seeds:
                for norder in norders:
                    for x0, y0, t0, beta in appsys:
                        waves_param.append({'dkx':dkx,
                                            'dky': dky,
                                            'nx':nx,
                                            'ny':ny,
                                            'ipol':ipol,
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
    dkx = request.param['dkx']
    dky = request.param['dky']
    nx = request.param['nx']
    ny = request.param['ny']
    ipol = request.param['ipol']
    seed = request.param['seed']
    norder = request.param['norder']
    x0 = request.param['x0']
    y0 = request.param['y0']
    t0 = request.param['t0']
    beta = request.param['beta']
    mysys = corsys.CorSys(x0, y0, t0, beta)
    cfuns, hfuns = tfun.create_Tfuns_2d(nx, ny, ipol, seed)
    swd_anal = shape_4.Shape4(dkx, dky, nx, ny, cfuns, hfuns, norder, mysys)
    swd_anal.tmpdir = str(tmp_path_factory.mktemp("swd"))
    file_swd = os.path.join(swd_anal.tmpdir, "shp_4.swd")
    swd_anal.write_swd(file_swd, dt=0.1, nsteps=11)
    # gfortran-9 does not allow a file to be concurrently opened twice...
    file_swd_copy = os.path.join(swd_anal.tmpdir, "shp_4_copy.swd")
    shutil.copy(file_swd, file_swd_copy)
    swd_nums = []
    for impl in [1, 2]:
        if impl == 1 or (impl == 2 and nx == ny and abs(dkx - dky) < 0.00001):
            # This trick of copy file should be removed when gfortran-10 is released
            if impl == 1:
                path = file_swd
            else:
                path = file_swd_copy
            swd_num = SpectralWaveData(path, x0, y0, t0, beta, rho=1025.0, impl=impl, dc_bias=True)
            swd_nums.append((impl, swd_num))

    yield swd_anal, swd_nums

    # Teardown code
    for impl, swd in swd_nums:
        swd.close()

xs = [0.0, 8.4]          # application x-positions for evaluation
ys = [0.0, -5.2]         # application y-positions for evaluation
zs = [0.0, -2.3, 1.7]    # application z-positions for evaluation
ts = [0.2, 0.23, 0.52]   # application time: First exact on 2dt, Then < dt, then > dt

def test_waves(make_waves):
    swd_anal, swd_nums = make_waves

    nx = swd_anal.nx
    ny = swd_anal.ny

    for jx in range(nx + 1):
        for jy in range(ny, ny + 1):
            swd_anal.dump_spectral_fun(jx, jy, dt=0.01, tmax=1.0)

    for t in ts:

        for impl, swd in swd_nums:
            swd.update_time(t)

        for x in xs:
            for y in ys:
                for z in zs:

                    phi_anal = swd_anal.phi(x, y, z, t)
                    for impl, swd in swd_nums:
                        phi_num = swd.phi(x, y, z)
                        assert math.isclose(phi_anal, phi_num, rel_tol=1e-04, abs_tol=1e-04)

                    stream_anal = swd_anal.stream(x, y, z, t)
                    for impl, swd in swd_nums:
                        stream_num = swd.stream(x, y, z)
                        assert math.isclose(stream_anal, stream_num, rel_tol=1e-04, abs_tol=1e-04)

                    phi_t_anal = swd_anal.phi_t(x, y, z, t)
                    for impl, swd in swd_nums:
                        phi_t_num = swd.phi_t(x, y, z)
                        assert math.isclose(phi_t_anal, phi_t_num, rel_tol=1e-04, abs_tol=1e-04)

                    grad_phi_anal = swd_anal.grad_phi(x, y, z, t)
                    for impl, swd in swd_nums:
                        grad_phi_num = swd.grad_phi(x, y, z)
                        assert math.isclose(grad_phi_num.x, grad_phi_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num.y, grad_phi_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_phi_num.z, grad_phi_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_anal = swd_anal.grad_phi_2nd(x, y, z, t)
                    for impl, swd in swd_nums:
                        g2nd_num = swd.grad_phi_2nd(x, y, z)
                        assert math.isclose(g2nd_num.xx, g2nd_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num.xy, g2nd_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num.xz, g2nd_anal['xz'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num.yy, g2nd_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num.yz, g2nd_anal['yz'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_num.zz, g2nd_anal['zz'], rel_tol=1e-04, abs_tol=1e-04)

                    acc_euler_anal = swd_anal.acc_euler(x, y, z, t)
                    for impl, swd in swd_nums:
                        acc_euler_num = swd.acc_euler(x, y, z)
                        assert math.isclose(acc_euler_num.x, acc_euler_anal['x'], rel_tol=1e-03, abs_tol=1e-03)
                        assert math.isclose(acc_euler_num.y, acc_euler_anal['y'], rel_tol=1e-03, abs_tol=1e-03)
                        assert math.isclose(acc_euler_num.z, acc_euler_anal['z'], rel_tol=1e-03, abs_tol=1e-03)

                    acc_particle_anal = swd_anal.acc_particle(x, y, z, t)
                    for impl, swd in swd_nums:
                        acc_particle_num = swd.acc_particle(x, y, z)
                        assert math.isclose(acc_particle_num.x, acc_particle_anal['x'], rel_tol=1e-03, abs_tol=1e-03)
                        assert math.isclose(acc_particle_num.y, acc_particle_anal['y'], rel_tol=1e-03, abs_tol=1e-03)
                        assert math.isclose(acc_particle_num.z, acc_particle_anal['z'], rel_tol=1e-03, abs_tol=1e-03)

                    prs_anal = swd_anal.pressure(x, y, z, t, rho=1025.0, grav=9.81)
                    for impl, swd in swd_nums:
                        prs_num = swd.pressure(x, y, z)
                        assert math.isclose(prs_num, prs_anal, rel_tol=1e-03)

                    elev_anal = swd_anal.elev(x, y, t)
                    for impl, swd in swd_nums:
                        elev_num = swd.elev(x, y)
                        assert math.isclose(elev_num, elev_anal, rel_tol=1e-04, abs_tol=1e-04)

                    elev_t_anal = swd_anal.elev_t(x, y, t)
                    for impl, swd in swd_nums:
                        elev_t_num = swd.elev_t(x, y)
                        assert math.isclose(elev_t_num, elev_t_anal, rel_tol=1e-04, abs_tol=1e-04)

                    grad_elev_anal = swd_anal.grad_elev(x, y, t)
                    for impl, swd in swd_nums:
                        grad_elev_num = swd.grad_elev(x, y)
                        assert math.isclose(grad_elev_num.x, grad_elev_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num.y, grad_elev_anal['y'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(grad_elev_num.z, grad_elev_anal['z'], rel_tol=1e-04, abs_tol=1e-04)

                    g2nd_elev_anal = swd_anal.grad_elev_2nd(x, y, t)
                    for impl, swd in swd_nums:
                        g2nd_elev_num = swd.grad_elev_2nd(x, y)
                        assert math.isclose(g2nd_elev_num.xx, g2nd_elev_anal['xx'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num.xy, g2nd_elev_anal['xy'], rel_tol=1e-04, abs_tol=1e-04)
                        assert math.isclose(g2nd_elev_num.yy, g2nd_elev_anal['yy'], rel_tol=1e-04, abs_tol=1e-04)

                    if abs(swd_anal.sys.t0 - 0.34) < 1.0e-6 and abs(t - 0.52) < 1.0e-6:
                        # Check if strip method produce a swd file with consistent kinematic results.
                        for impl, swd in swd_nums:
                            file_swd_strip = os.path.join(swd_anal.tmpdir, 'stripped.swd')
                            swd.strip(tmin=0.31, tmax=0.61, file_swd=file_swd_strip)
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
                    for impl, swd in swd_nums:
                        swd.convergence(x, y, z, csv=file_convergence)
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

        depth_anal = swd_anal.bathymetry(x_app=3.4, y_app=-5.3)
        nvec_floor_anal = swd_anal.bathymetry_nvec(x_app=3.4, y_app=-5.3)
        for impl, swd in swd_nums:
            depth_num = swd.bathymetry(x=3.4, y=-5.3)
            if depth_anal < 0:
                assert depth_num < 0
            else:
                assert math.isclose(depth_num, depth_anal, rel_tol=1e-5)
            nvec_floor_num = swd.bathymetry_nvec(x=3.4, y=-5.3)
            assert math.isclose(nvec_floor_anal['x'], nvec_floor_num.x, abs_tol=1e-5)
            assert math.isclose(nvec_floor_anal['y'], nvec_floor_num.y, abs_tol=1e-5)
            assert math.isclose(nvec_floor_anal['z'], nvec_floor_num.z, rel_tol=1e-5)

        for impl, swd in swd_nums:
            # Check the get method for int, float and strings arguments
            nx_file = swd.get('nx')
            assert isinstance(nx_file, int)
            assert nx_file ==  swd_anal.nx

            dt_swd = swd.get('dt')
            assert isinstance(dt_swd, float)
            assert math.isclose(dt_swd, 0.1, rel_tol=1e-5)

            cid = swd.get('cid')
            assert isinstance(cid, str)
            assert cid == "{'ole':1, 'dole':2, 'doffen':3}"
