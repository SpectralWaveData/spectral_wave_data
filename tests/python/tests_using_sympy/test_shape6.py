"""
Test spectra_wave_shape_2 against analytical fields evaluated by sympy
"""

import sys, os
import math, cmath

import pytest
import numpy as np

import shape_6
import corsys
import tfun

from spectral_wave_data import SpectralWaveData
from spectral_wave_data.tools import airy

assert sys.version_info > (3, 4)

# We check all permutations of...

fields = []
fields.append({'grav': 9.81,
               'depth': -1,
               'amps': [2.0, 0.0],
               'dirs': [0.0, 0.0],
               'phases': [0.0, 90.0],
               'Twaves': {5.3, 13.2}
             })

fields.append({'grav': 9.81,
               'depth': 6.0,
               'amps': [2.0, 3.0],
               'dirs': [122.0, 62.0],
               'phases': [0.0, 90.0],
               'Twaves': {5.3, 13.2}
             })

fields.append({'grav': 9.81,
               'depth': 6.0,
               'amps': [2.0, 3.0],
               'dirs': [0.0, 0.0],
               'phases': [24.0, -120.0],
               'Twaves': {5.3, 13.2}
             })


#norders = [-1, 3]   # Taylor expansion order in z-direction
norders = [-1, 0, 1, 2]       # Different schemes for handling z>0
appsys = [(0.0, 0.0, 0.0, 0.0),  #'x0, y0, t0, beta'
          (13.8, 0.0, 0.0, 0.0),
          (0.0, 18.2, 0.0, 0.0),
          (0.0, 0.0, 0.34, 0.0),
          (-5.3, 18.2, 0.0, 37.4)] # Coordinate transformations

waves_param = []
for id_field in range(len(fields)):
    for norder in norders:
        for x0, y0, t0, beta in appsys:
            waves_param.append({'field':id_field,
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
    field = fields[request.param['field']]
    grav = field['grav']
    depth = field['depth']
    amps = np.array(field['amps'])
    dirs = np.array(field['dirs']) * np.pi / 180.0
    phases = np.array(field['phases']) * np.pi / 180.0
    Twaves = field['Twaves']
    kwaves = [airy.omega2kwave((2*np.pi/T), depth, grav) for T in Twaves]
    norder = request.param['norder']
    x0 = request.param['x0']
    y0 = request.param['y0']
    t0 = request.param['t0']
    beta = request.param['beta']
    mysys = corsys.CorSys(x0, y0, t0, beta)
    swd_anal = shape_6.Shape6(amps, dirs, phases, kwaves, depth, norder, grav, mysys)
    swd_anal.tmpdir = str(tmp_path_factory.mktemp("swd"))
    file_swd = os.path.join(swd_anal.tmpdir, "shp_6.swd")
    swd_anal.write_swd(file_swd)
    swd_num = SpectralWaveData(file_swd, x0, y0, t0, beta, rho=1025.0, impl=0, norder=norder, dc_bias=True)

    yield swd_anal, swd_num

    # Teardown code
    swd_num.close()

xs = [0.0, 8.4]          # application x-positions for evaluation
ys = [0.0, -5.2]         # application y-positions for evaluation
zs = [0.0, -2.3, 3.7]    # application z-positions for evaluation
ts = [0.0, 15.2]         # application time

def test_waves(make_waves):
    swd_anal, swd_num = make_waves

    for t in ts:

        swd_num.update_time(t)

        for x in xs:
            for y in ys:
                for z in zs:

                    phi_anal = swd_anal.phi(x, y, z, t)
                    phi_num = swd_num.phi(x, y, z)
                    assert math.isclose(phi_anal, phi_num, rel_tol=1e-04, abs_tol=1e-04)

                    stream_anal = swd_anal.stream(x, y, z, t)
                    stream_num = swd_num.stream(x, y, z)
                    assert math.isclose(stream_anal, stream_num, rel_tol=1e-04, abs_tol=1e-04)

                    phi_t_anal = swd_anal.phi_t(x, y, z, t)
                    phi_t_num = swd_num.phi_t(x, y, z)
                    assert math.isclose(phi_t_anal, phi_t_num, rel_tol=1e-04, abs_tol=1e-04)

                    grad_phi_anal = swd_anal.grad_phi(x, y, z, t)
                    grad_phi_num = swd_num.grad_phi(x, y, z)
                    assert math.isclose(grad_phi_num.x, grad_phi_anal['x'], rel_tol=1e-04, abs_tol=1e-04)
                    assert math.isclose(grad_phi_num.y, grad_phi_anal['y'], rel_tol=1e-04, abs_tol=1e-04)

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

        depth_anal = swd_anal.bathymetry(x_app=3.4, y_app=-5.3)
        nvec_floor_anal = swd_anal.bathymetry_nvec(x_app=3.4, y_app=-5.3)
        depth_num = swd_num.bathymetry(x=3.4, y=-5.3)
        if depth_anal < 0:
            assert depth_num < 0
        else:
            assert math.isclose(depth_num, depth_anal, rel_tol=1e-5)
        nvec_floor_num = swd_num.bathymetry_nvec(x=3.4, y=-5.3)
        assert math.isclose(nvec_floor_anal['x'], nvec_floor_num.x, abs_tol=1e-5)
        assert math.isclose(nvec_floor_anal['y'], nvec_floor_num.y, abs_tol=1e-5)
        assert math.isclose(nvec_floor_anal['z'], nvec_floor_num.z, rel_tol=1e-5)

        # Check the get method for int, float and strings arguments
        n = swd_num.get('n')
        assert isinstance(n, int)
        assert n ==  swd_anal.n

        grav = swd_num.get('grav')
        assert isinstance(grav, float)
        assert math.isclose(grav, swd_anal.grav, rel_tol=1e-5)

        swd_class = swd_num.get('class')
        assert isinstance(swd_class, str)
        assert swd_class == "spectral_wave_data_shape_6_impl_1"
