# Build test_interpolation programs

cmake_minimum_required(VERSION 3.4.1)

project(test_interpolation Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(DIR_API "${CMAKE_SOURCE_DIR}/../../../../src/api/fortran")

set(SRC_API
  ${DIR_API}/kind_values.f90
  ${DIR_API}/spectral_interpolation.f90)
source_group(API FILES ${SRC_API})

set(SRC_TEST0
  ../test_spectral_interpolation_ipol_0.f90)

set(SRC_TEST1
  ../test_spectral_interpolation_ipol_1.f90)

set(SRC_TEST2
  ../test_spectral_interpolation_ipol_2.f90)

add_executable(test0 ${SRC_API} ${SRC_TEST0})
add_executable(test1 ${SRC_API} ${SRC_TEST1})
add_executable(test2 ${SRC_API} ${SRC_TEST2})

# Set Visual Studio StartUp Project
set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT test0)
