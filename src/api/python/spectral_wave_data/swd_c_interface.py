# -*- coding: utf-8 -*-
    
"""
:platform: Linux, Windows, python 2.7 and 3.x
:synopsis: Defines the Python-C interface of spectral_wave_data

Author  - Jens Bloch Helmers, DNVGL
Created - 2019-08-11
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import platform
from ctypes import c_bool, c_double, c_int, c_char_p, c_void_p, Structure, CDLL, POINTER
from .ISO_Fortran_binding import CFI_cdesc_t_1D, CFI_cdesc_t_2D, CFI_cdesc_t_3D

assert sys.version_info >= (2, 7, 11)

HERE = os.path.dirname(os.path.abspath(__file__))
if platform.system() == 'Windows':
    # The Intel Redistributable Libraries are expected to be statical linked to the dll.
    # Consequently, we skip this:
    # intel_redist_path = os.getenv('INTEL_DEV_REDIST')
    # if intel_redist_path is None:
    #     msg = """
    #           To apply "spectral_wave_data" you need to install the latest version of:
    #           1) Redistributable Libraries for Intel® C++ and Fortran Compilers for Windows
    #              This package can be downloaded for free from intel.com
    #           2) You also need Microsoft Visual C++ redistributable or
    #              Visual Studio with C++ tools and libraries.
    #              These tools can be downloaded from microsoft.com
    #           """
    #     raise AssertionError(msg)
    # if sys.version_info >= (3, 8):
    #     intel_redist_path = os.path.join(intel_redist_path, 'redist', 'intel64', 'compiler')
    #     os.add_dll_directory(HERE)
    #     os.add_dll_directory(intel_redist_path)
    swdlib = CDLL(str(os.path.join(HERE, 'SpectralWaveData.dll')))
elif platform.system() == 'Linux':
    swdlib = CDLL(str(os.path.join(HERE, 'libSpectralWaveData.so')))
else:
    raise AssertionError('Not supported platform: ' + platform.system())

"""
================================================================================================
BEGIN interface definition to the C-implementation
================================================================================================
NOTE: STRANGE ERRORS may occur if this interface does not comply with the original C source code.
"""

class vecswd(Structure):
     _fields_ = [("x", c_double), ("y", c_double), ("z", c_double)]


class vecphi2ndswd(Structure):
    _fields_ = [("xx", c_double), ("xy", c_double), ("xz", c_double),
                ("yy", c_double), ("yz", c_double), ("zz", c_double)]

class vecelev2ndswd(Structure):
    _fields_ = [("xx", c_double), ("xy", c_double), ("yy", c_double)]

swdlib.swd_api_allocate.argtypes = [c_char_p, c_double, c_double,
                                    c_double, c_double, c_double,
                                    c_int, c_int, c_int, c_int,
                                    c_int, c_bool]
swdlib.swd_api_allocate.restype = c_void_p

swdlib.swd_api_update_time.argtypes = [c_void_p, c_double]
swdlib.swd_api_update_time.restype = c_void_p

swdlib.swd_api_phi.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_phi.restype = c_double

swdlib.swd_api_stream.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_stream.restype = c_double

swdlib.swd_api_phi_t.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_phi_t.restype = c_double

swdlib.swd_api_grad_phi.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_grad_phi.restype = vecswd

swdlib.swd_api_grad_phi_2nd.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_grad_phi_2nd.restype = vecphi2ndswd

swdlib.swd_api_acc_euler.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_acc_euler.restype = vecswd

swdlib.swd_api_acc_particle.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_acc_particle.restype = vecswd

swdlib.swd_api_elev.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_elev.restype = c_double

swdlib.swd_api_elev_t.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_elev_t.restype = c_double

swdlib.swd_api_grad_elev.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_grad_elev.restype = vecswd

swdlib.swd_api_grad_elev_2nd.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_grad_elev_2nd.restype = vecelev2ndswd

swdlib.swd_api_pressure.argtypes = [c_void_p, c_double, c_double, c_double]
swdlib.swd_api_pressure.restype = c_double

swdlib.swd_api_bathymetry.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_bathymetry.restype = c_double

swdlib.swd_api_bathymetry_nvec.argtypes = [c_void_p, c_double, c_double]
swdlib.swd_api_bathymetry_nvec.restype = vecswd

swdlib.swd_api_convergence.argtypes = [c_void_p, c_double, c_double, c_double, c_char_p]
swdlib.swd_api_convergence.restype = c_void_p

swdlib.swd_api_strip.argtypes = [c_void_p, c_double, c_double, c_char_p]
swdlib.swd_api_strip.restype = c_void_p

swdlib.swd_api_get_chr.argtypes = [c_void_p, c_char_p]
swdlib.swd_api_get_chr.restype = c_char_p

swdlib.swd_api_get_int.argtypes = [c_void_p, c_char_p]
swdlib.swd_api_get_int.restype = c_int

swdlib.swd_api_get_bool.argtypes = [c_void_p, c_char_p]
swdlib.swd_api_get_bool.restype = c_bool

swdlib.swd_api_get_real.argtypes = [c_void_p, c_char_p]
swdlib.swd_api_get_real.restype = c_double

swdlib.swd_api_error_raised.argtypes = [c_void_p]
swdlib.swd_api_error_raised.restype = c_bool

swdlib.swd_api_error_get_id.argtypes = [c_void_p]
swdlib.swd_api_error_get_id.restype = c_int

swdlib.swd_api_error_get_msg.argtypes = [c_void_p]
swdlib.swd_api_error_get_msg.restype = c_char_p

swdlib.swd_api_error_clear.argtypes = [c_void_p]
swdlib.swd_api_error_clear.restype = c_void_p

swdlib.swd_api_close.argtypes = [c_void_p]
swdlib.swd_api_close.restype = c_void_p

swdlib.swd_api_elev_fft.argtypes = [c_void_p, c_int, c_int]
swdlib.swd_api_elev_fft.restype = POINTER(CFI_cdesc_t_2D)

swdlib.swd_api_grad_phi_fft.argtypes = [c_void_p, c_double, c_int, c_int]
swdlib.swd_api_grad_phi_fft.restype = POINTER(CFI_cdesc_t_3D)

swdlib.swd_api_x_fft.argtypes = [c_void_p, c_int]
swdlib.swd_api_x_fft.restype = POINTER(CFI_cdesc_t_1D)

swdlib.swd_api_y_fft.argtypes = [c_void_p, c_int]
swdlib.swd_api_y_fft.restype = POINTER(CFI_cdesc_t_1D)

"""
================================================================================================
END interface to C implementation
================================================================================================
"""
