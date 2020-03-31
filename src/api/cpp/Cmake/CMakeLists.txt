# Build a simple example of C++-program using
# the Fortran spectral-wave-data static library which is build first

cmake_minimum_required(VERSION 3.4.1)

set(PROG swd_cpp_app_example)
set(LIB spectral_wave_data)

project(${PROG} CXX Fortran)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_API_C ../../c)
set(DIR_SRC_API_CPP ..)
set(DIR_SRC_APP ../../../api_examples)

set(SRC_API_F
    ${DIR_SRC_API_F}/open_swd_file.f90
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
    ${DIR_SRC_API_CPP}/kind_values.f90     # Note we apply from the C++ directory
)

set(SRC_CPP
    ${DIR_SRC_APP}/application_swd_api.cpp
    ${DIR_SRC_API_CPP}/SpectralWaveData.cpp
    ${DIR_SRC_API_CPP}/SpectralWaveData.h
    ${DIR_SRC_API_C}/spectral_wave_data.h
)

add_library(${LIB} STATIC ${SRC_API_F})
target_include_directories(${LIB} INTERFACE 
                           ${DIR_SRC_API_C}
                           ${DIR_SRC_API_CPP})

add_executable(${PROG} ${SRC_CPP})
target_link_libraries(${PROG} PRIVATE ${LIB}) 
