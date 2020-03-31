# -*- coding: utf-8 -*-
    
"""
THIS IS A MOCK OBJECT TO PREVENT IMPORT OF C-FUNCTIONS
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import sys
assert sys.version_info > (2, 6)

class vecswd(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class vecphi2ndswd(object):
    def __init__(self, xx, xy, xz, yy, yz, zz):
        self.xx = xx
        self.xy = xy
        self.xz = xz
        self.yy = yy
        self.yz = yz
        self.zz = zz

class vecelev2ndswd(object):
    def __init__(self, xx, xy, yy):
        self.xx = xx
        self.xy = xy
        self.yy = yy

class Cswd(object):
    def __init__(self):
        self.dat = None

def swd_api_allocate(swd, file_swd, x0, y0, t0, beta, rho, nsumx, nsumy,
                     impl, ipol, norder, dc_bias):
    return Cswd()

def swd_api_update_time(swd, time):
    return None

def swd_api_phi(swd, x, y, z):
    return 0.0

def swd_api_stream(swd, x, y, z):
    return 0.0

def swd_api_phi_t(swd, x, y, z):
    return 0.0

def swd_api_grad_phi(swd, x, y, z):
    res = vecswd(x, y, z)
    return res

def swd_api_grad_phi_2nd(swd, x, y, z):
    res = vecphi2ndswd(x, y, z, x, y, z)
    return res

def swd_api_acc_euler(swd, x, y, z):
    res = vecswd(x, y, z)
    return res

def swd_api_acc_particle(swd, x, y, z):
    res = vecswd(x, y, z)
    return res

def swd_api_elev(swd, x, y):
    return 0.0

def swd_api_elev_t(swd, x, y):
    return 0.0

def swd_api_grad_elev(swd, x, y):
    res = vecswd(x, y, z)
    return res

def swd_api_grad_elev_2nd(swd, x, y):
    res = vecelev2ndswd(x, y, z)
    return res

def swd_api_pressure(swd, x, y, z):
    return 0.0

def swd_api_bathymetry(swd, x, y):
    return 0.0

def swd_api_bathymetry_nvec(swd, x, y):
    res = vecswd(x, y, z)
    return res

def swd_api_convergence(swd, x, y, z, csv):
    return None

def swd_api_strip(swd, tmin, tmax, file_swd):
    return None

def swd_api_get_chr(swd, name):
    return b"nothing"

def swd_api_get_int(swd, name):
    return 10

def swd_api_get_bool(swd, name):
    return False

def swd_api_get_real(swd, name):
    return 10.0

def swd_api_error_raised(swd):
    return False

def swd_api_error_get_id(swd):
    return 0

def swd_api_error_get_msg(swd):
    return b"something"

def swd_api_error_clear(swd):
    return None

def swd_api_close(swd):
    return None
