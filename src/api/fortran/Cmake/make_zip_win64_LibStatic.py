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

build_dir = "Build_Win64_LibStatic"
release_dir = os.path.join(build_dir, 'Release')

zip_name = 'swd_lib_static_f_win64_%s' % version_full
zip_dir = os.path.join(build_dir, zip_name)
lib_dir = os.path.join(zip_dir, 'lib')
mod_dir = os.path.join(zip_dir, 'mod')

if os.path.isdir(zip_dir):
    shutil.rmtree(zip_dir)
os.mkdir(zip_dir)

os.mkdir(lib_dir)
for path in glob.glob(os.path.join(release_dir, '*.lib')):
    shutil.copy(path, lib_dir)

os.mkdir(mod_dir)
for path in glob.glob(os.path.join(release_dir, '*.mod')):
    shutil.copy(path, mod_dir)

readme = """
This is a binary distribution as part the spectral_wave_data library.
---------------------------------------------------------------------

spectral_wave_data version = %s

Requirements:
=============

  1) To be applied in a Fortran-program.
  2) Your operating system is Windows-10
  3) You apply x64 architecture for your compilation.
  4) Your Fortran-compiler is binary compatible with recent versions
     of Intel Fortran.
  5) Your Fortran-program applies double precision in all function 
     signatures interacting with this library. 

The 'lib' folder contains the static library of Fortran code.

The 'mod' folder contains the Fortran module interfaces.

This distribution is covered by the copyright and license
as described in the GitHub repository: 
https://github.com/SpectralWaveData/spectral_wave_data.
""" % (version_full)

with open(os.path.join(zip_dir, 'README.txt'), 'w') as f:
    f.write(readme)

shutil.copy(os.path.join(swd_root, 'LICENSE.txt'), zip_dir)

shutil.make_archive(zip_name, 'zip', zip_dir)

print("Your archive is created:")
print(zip_name + '.zip')
