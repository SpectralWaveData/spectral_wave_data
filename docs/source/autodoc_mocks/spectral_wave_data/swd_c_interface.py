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

import spectral_wave_data.swdlib as swdlib
from .swdlib import vecswd, vecphi2ndswd, vecelev2ndswd
