# -*- coding: utf-8 -*-

"""
:platform: Linux, Windows, python 3.x
:synopsis: Define version number of spectral-wave-data.

Author  - Jens Bloch Helmers, DNV
Created - 2020-01-24
"""

from __future__ import print_function
from __future__ import unicode_literals

import os
import re

#=====================================================================================================
# TO UPDATE THE SWD SOURCE CODE YOU MUST RUN THIS PYTHON SCRIPT AFTER YOU MODIFY THE TAGS:
tag_major = 1          # int  (Only increment if incompatibilities are introduced)
tag_minor = 0          # int  (Increment if new features are introduced. 100% backward compatibility)
tag_patch = 0          # int  (Increment if bug fix, improved documentation etc. No new features)
tag_prerel = 'rc4'     # E.g. 'beta.1', 'beta.1.debug.10', 'rc1'
tag_build = ''         # E.g. 'VendorX', branch and/or signature. (dot as separation)
# No need to edit below unless you really have to do it...
#=====================================================================================================

# Note that both the Sphinx-documentation and the Python-whl-package setup get
# the relevant version number by importing this Python module directly. For compiled languages
# like Fortran this script creates relevant source code containing the new version number.

# Construct standard conforming SEMVER (https://semver.org/) version numbers to be imported in software
version_short = '%i.%i' % (tag_major, tag_minor)
version_core = '%i.%i.%i' % (tag_major, tag_minor, tag_patch)
version_full = version_core
if tag_prerel != '':
    version_full += '-' + tag_prerel
if tag_build != '':
    version_full += '+' + tag_build

# Check strict SEMVER-rules for legal version number (The compile string is copied from semver.org.)
p = re.compile(r'^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)'
               '(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)'
               '(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?'
               '(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$')
assert p.match(version_core)
assert p.match(version_full)

def make_fortran_version():
    """Implicit also C, C++ and Python"""
    path_fortran = os.path.join('src', 'api', 'fortran', 'swd_version.f90')
    out = open(path_fortran, 'w')
    out.write("module swd_version\n\n")
    out.write("! Do not edit this file!!!\n")
    out.write("! This file is made by running 'version_definition.py' in the root directory.\n\n")
    out.write("character(len=*), parameter :: version = '%s'\n\n" % version_full)
    out.write("end module swd_version\n")
    out.close()
    print("A new file '%s' is created. You should recompile your SWD libraries..." % path_fortran)

def make_python_version():
    """Same version as C, C++ and Fortran"""
    path = os.path.join("src", "api", "python", "distribution", "VERSION.dat")
    out = open(path, 'w')
    out.write("%s\n" % version_full)
    out.close()

if __name__ == "__main__":
    make_fortran_version()
    make_python_version()
