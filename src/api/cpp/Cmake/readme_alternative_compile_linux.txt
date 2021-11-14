#!/bin/bash

# Some prefer this recipe on Linux (see comments below):
rmdir /s /q Build_Linux
mkdir Build_Linux
cp CMakeLists.txt Build_Linux
cd Build_Linux
export FC=/usr/bin/gfortran-10
export CC=/usr/bin/gcc-10
export CXX=/usr/bin/g++-10
export CPP=/usr/bin/cpp-10
cmake -G "Unix Makefiles" ..
make


# Requires gfortran-10:   apt-get install gfortran-10 (May still work on gcc-9, but not gcc-8)
# Should perhaps add the line
# set(CMAKE_POSITION_INDEPENDENT_CODE ON) 
# in CMakeLists.txt
#
# Instructions:
# Put above commands in a script (e.g. make_linux.sh) in 
# spectral-wave-data/src/api/c/Cmake and 
# run it by ./make_linux.sh, making sure that 
# 'Allow executing file as program' is ticked off in permissions (Cmake GUI).
# This will make the libspectral_wave_data.a static library.
# To use it in a C-program, i.e. put the libspectral_wave_data.a 
# file, the spectral_wave_data.h file and the program (.c) file
# in the same directory.
# Add the two lines
# #include <stdbool.h>
# #include "spectral_wave_data.h"    // Namespace convention: swd_api_*
# at the top of the .c file
# and compile by
# gcc my_program.c -o my_program -L/usr/lib/gcc/x86_64-linux-gnu/10 libspectral_wave_data.a -lgfortran -lm -lquadmath -lm -O2 
