## Software Prerequisites

### Windows-10 or 11

On Windows the precompiled libraries and Python install wheels are build
using Intel ifort OneAPI 2023.1 and Microsoft Visual C++ 2022. 

If you encounter error messages related to loading the Python package,
you have a corrupt setup of your Python environment. The best solution is 
to install spectral_wave_data in a new Python environment.

### Linux

There are many different Linux distributions using different incompatible
versions of system libraries. The current python distribution on PyPI is 
compiled using Intel OneAPI to avoid incompatible requirements of 
gcc/g++/gfortran libraries. (PyPI accepts only very old versions of GNU
libraries. No modern versions are accepted.)

The uploaded Linux binaries (not Python wheels) are compiled and tested 
using gfortran/gcc 11.x on the Ubuntu x86_64 architecture.
They are therefore expected to work on related distroes. However, in case
of problems when loading system libraries, you need to install the 
Gnu Compiler Collection (gcc), including gfortran, version 11.x or newer.

If you have another Linux distro you may still try the precompiled libraries.
If they don't work, it is not difficult to compile from source:

- Install and apply gcc/g++/gfortran version 11.x or newer
- Run the [Cmake](https://cmake.org/) scripts included in this repository.
