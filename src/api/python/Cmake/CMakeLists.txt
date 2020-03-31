# Build SpectralWaveData library for Python package

cmake_minimum_required(VERSION 3.4.1)

set(LIB SpectralWaveData)

project(${LIB} C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_BUILD_TYPE Release)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_API_C ../../c)
set(DIR_SRC_API_PY ..)

set(SRC_API_F
    ${DIR_SRC_API_F}/open_swd_file.F90
    ${DIR_SRC_API_F}/spectral_interpolation.f90
    ${DIR_SRC_API_F}/spectral_wave_data.f90
    ${DIR_SRC_API_F}/spectral_wave_data_allocate.f90
    ${DIR_SRC_API_F}/spectral_wave_data_c.f90
    ${DIR_SRC_API_F}/spectral_wave_data_error.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_1_impl_1.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_2_impl_1.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_3_impl_1.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_4_impl_1.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_4_impl_2.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_5_impl_1.f90
    ${DIR_SRC_API_F}/spectral_wave_data_shape_6_impl_1.f90
    ${DIR_SRC_API_F}/swd_version.f90
    ${DIR_SRC_API_PY}/kind_values.f90     # Note we apply from the Python-directory
)

add_library(${LIB} SHARED ${SRC_API_F})
target_include_directories(${LIB} INTERFACE ${DIR_SRC_API_C})

