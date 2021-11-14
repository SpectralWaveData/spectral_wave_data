This Python package provides the spectral_wave_data ocean wave model.

Programmed by: Jens B. Helmers & Odin Gramstad, DNVGL

Installation:
------------------------------------------------------------------

0) Ensure the following is already installed on your computer:
    a) Python 2.7 or 3.x
    b) Python package wheel installed (recent version)

1) From the command prompt run:
   pip install spectral_wave_data-xxx-correct-platform.whl

   The spectral_wave_data package can be uninstalled using the command:
   pip uninstall spectral_wave_data

   To check actual version of library
   pip show spectral_wave_data

------------------------------------------------------------------


Notes for building SWD(Python) for windows:
===========================================

To avoid end users dependency on Intel Redistributable
https://software.intel.com/en-us/articles/intel-compilers-redistributable-libraries-by-version

it is recommended to build the Fortran library SpectralWaveData.dll using the /MT
Intel compiler flag. Then all required Intel routines will be included in the above dll
and the Python whl file.


Notes for building SWD(Python) for Linux:
=========================================

You may need libraries included by gfortran-10/gcc-10 or later. 
(gfortran/9 and gcc-9 may still work but not version 8)

If you want to create a manylinux wheel, the "easiest" is to get an old
version of linux, install Intel Fortran and compile the libSpectralWaveData.so
file with ifort and the "-static-intel" flag to make it depend only on the
(old) system files. You can then run the "auditwheel" utility to check whl
file compatibility and make the wheel file "universal" (manylinux-compatible).
Getting gfortran version 9+ to run on the old manylinux images is very hard.
