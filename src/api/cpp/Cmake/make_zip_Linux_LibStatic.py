import os
import sys
import glob
import shutil

"""
Collect compiled components and store them in a zip-file
to be uploaded as a part of the GitHub binary release.
"""

# To read version number from version_definition.py in root directory
swd_root = os.path.abspath(os.path.join('..', '..', '..', '..'))
sys.path.insert(0, swd_root)
from version_definition import version_full

build_dir = "Build_Linux_LibStatic"

zip_name = 'swd_lib_static_cpp_linux_%s' % version_full
zip_dir = os.path.join(build_dir, zip_name)
lib_dir = os.path.join(zip_dir, 'lib')
mod_dir = os.path.join(zip_dir, 'mod')
inc_dir = os.path.join(zip_dir, 'inc')
cpp_dir = os.path.join(zip_dir, 'cpp')
src_c_dir =  os.path.join('..', '..', 'c')
src_cpp_dir =  os.path.join('..')

if os.path.isdir(zip_dir):
    shutil.rmtree(zip_dir)
os.mkdir(zip_dir)

os.mkdir(lib_dir)
for path in glob.glob(os.path.join(build_dir, '*.a')):
    shutil.copy(path, lib_dir)

os.mkdir(mod_dir)
for path in glob.glob(os.path.join(build_dir, '*.mod')):
    shutil.copy(path, mod_dir)

os.mkdir(inc_dir)
shutil.copy(os.path.join(src_c_dir, 'spectral_wave_data.h'), inc_dir)
shutil.copy(os.path.join(src_cpp_dir, 'SpectralWaveData.h'), inc_dir)

os.mkdir(cpp_dir)
shutil.copy(os.path.join(src_cpp_dir, 'SpectralWaveData.cpp'), cpp_dir)


readme = """
This is a binary distribution as part the spectral_wave_data library.
---------------------------------------------------------------------

spectral_wave_data version = %s

Requirements:
=============

  1) To be applied in a C++ program.
  2) Your operating system is Linux
  3) You apply x64 architecture for your compilation.
  4) Your C++ compiler is binary compatible with recent versions
     of gcc/g++.
  5) Your C++ program applies double (not float) in all function 
     signatures interacting with this library. 

The 'lib' folder contains the static library of Fortran code.

The 'mod' folder contains the Fortran module interfaces.

The 'inc' folder contains the relevant C and C++ header files to include
in your C++ program.

The 'cpp' folder contains C++ code to be included in your C++ program.

This distribution is covered by the copyright and license
as described in the GitHub repository: 
https://github.com/SpectralWaveData/spectral_wave_data.
""" % (version_full)

with open(os.path.join(zip_dir, 'README.txt'), 'w') as f:
    f.write(readme)

shutil.copy(os.path.join(swd_root, 'LICENSE.txt'), zip_dir)

shutil.make_archive(zip_name, 'gztar', zip_dir)

print("Your archive is created:")
print(zip_name + '.tar.gz')
