This Python package provides the spectral_wave_data ocean wave model.

Programmed by: Jens B. Helmers & Odin Gramstad, DNVGL

Installation:
------------------------------------------------------------------

0) Ensure the following is already installed on your computer:
    a) Python 2.7 or 3.x
    b) Python package wheel installed

1) From the command prompt run:
   pip install spectral_wave_data-xxx-correct-plattform.whl

   The spectral_wave_data package can be unistalled using the command:
   pip uninstall spectral_wave_data

   To check actual version of library
   pip show spectral_wave_data

------------------------------------------------------------------


Notes for windows users:
========================

Unless you have installed a recent Intel Fortran and Microsoft C++ compiler you
may have to install some redistributable packages to be downloaded from Intel and Microsoft:

https://software.intel.com/en-us/articles/intel-compilers-redistributable-libraries-by-version
https://www.microsoft.com/en-us/download/details.aspx?id=52685


Notes for Linux users:
======================

You may need libraries included by gfortran-9/gcc-9 or later.

If you want to create a manylinux wheel, the "easiest" is to get an old
version of linux, install Intel Fortran and compile the libSpectralWaveData.so
file with ifort and the "-static-intel" flag to make it depend only on the
(old) system files. You can then run the "auditwheel" utility to check whl
file compatibility and make the wheel file "universal" (manylinux-compatible).
Getting gfortran version 9+ to run on the old manylinux images is very hard.
