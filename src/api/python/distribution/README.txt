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

You may need libraries included by gfortran-11/gcc-11 or later. 
(older versions may still work but not version 8 and below)

The new Intel OneAPI (free) should work fine. (Applied in the PyPI distribution of swd)
