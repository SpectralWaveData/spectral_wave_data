## Software Prerequisites

### Windows-10

On Windows the precompiled libraries and Python install wheels are build
using Intel Fortran 2020. Consequently, our software requires recent
core system DLLs provided by Intel. Hence, you should download and install the
[Intel redistributable library for Fortran 2020](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2020-compilers-for-windows)
using default install options.

If you already have a recent Intel Fortran compiler installed, you may skip
this installation. However, if you encounter some DLL loading problems, you 
must install the actual redistributable library. 

The above instructions relate to application of the Python distribution too.

If you still encounter error messages related to loading the Python module,
you have a corrupt setup of your general Python system. This could be identified
by the order of loading python packages. 
The best solution is to reinstall your Python system in a proper way.
However, a quick-fix is usually to load spectral_wave_data before any other
packages because spectral_wave_data requires a more recent version of some 
Intel DLLs on Windows. If some other package loads an older system DLL first,
you are lost. This rare issue is not present if you apply recommended 
procedures for 
setting up your Python system. Hence, consistently using e.g. pyenv or 
Anaconda environments. In general, manual fumbling of system paths is 
always a source of subtle problems and should be avoided.


### Linux

There are many different Linux distributions using different incompatible
versions of system libraries. Without using very old versions of gcc/gfortran
it is not possible to provide binaries working on all Linux distributions.

Since spectral_wave_data **requires features from version 9.x or newer** of 
the Gnu Compiler Collection (gcc), it is not possible to create a unified 
binary. This is why you don't find a Linux version of spectral_wave_data on
[PyPI](https://pypi.org/project/spectral-wave-data/).

The uploaded Linux binaries and Python wheels are compiled and tested 
using gfortran/gcc 9.x on the Ubuntu x86_64 architecture.
It is therefore expected to work on related distroes.

These binaries have also briefly been tested on CentOS-8 with no errors in
the pytest package. Consequently, they are also expected to work on 
RedHat and Fedora related distros.

If you on above systems encounter problems of loading some system 
libraries, you need to install the Gnu Compiler Collection
(gcc), including gfortran, version 9.x or newer on your system.

If you have another Linux distro you may still try the precompiled libraries.
If they don't work, it is not difficult to compile from source:

- Install and apply gcc/g++/gfortran version 9.x or newer
- Run the [Cmake](https://cmake.org/) scripts included in this repository.
- For Python, finally run the included bash scripts to compile a wheel
  install file for your specific system.

