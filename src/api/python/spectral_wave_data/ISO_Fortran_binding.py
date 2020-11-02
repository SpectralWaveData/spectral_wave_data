# -*- coding: utf-8 -*-
    
"""
:platform: Linux, Windows, python 2.7 and 3.x
:synopsis: The type definitions corresponding to ISO_Fortran_binding.h

Author  - Odin Gramstad, DNV GL
Created - 2020-10-30
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import platform, sys, os
from ctypes import (c_bool, c_double, c_char_p, c_void_p, 
                    c_int, c_int8, c_int16, c_int32, c_int64, c_size_t,
                    Structure, cast, POINTER, sizeof)

assert sys.version_info >= (2, 7, 11)

c_int8_t = c_int8
c_int16_t = c_int16
if sizeof(c_void_p) == 4:
     c_ptrdiff_t = c_int32
elif sizeof(c_void_p) == 8:
     c_ptrdiff_t = c_int64

# CFI type definitions
CFI_index_t = c_ptrdiff_t
CFI_rank_t = c_int8_t
CFI_attribute_t = c_int8_t
CFI_type_t = c_int16_t

# CFI_dim_t structure
class CFI_dim_t(Structure):
     _fields_ = [("lower_bound", CFI_index_t), 
                 ("extent", CFI_index_t), 
                 ("sm", CFI_index_t)]
                 
# CFI_cdesc_t structure for 1D array
class CFI_cdesc_t_1D(Structure):
    _fields_ = [("base_addr", c_void_p), 
                ("elem_len", c_size_t), 
                ("version", c_int),
                ("rank", CFI_rank_t),
                ("attribute", CFI_attribute_t),
                ("type", CFI_type_t),
                ("dim", CFI_dim_t*1)]

# CFI_cdesc_t structure for 2D array
class CFI_cdesc_t_2D(Structure):
    _fields_ = [("base_addr", c_void_p), 
                ("elem_len", c_size_t), 
                ("version", c_int),
                ("rank", CFI_rank_t),
                ("attribute", CFI_attribute_t),
                ("type", CFI_type_t),
                ("dim", CFI_dim_t*2)]
