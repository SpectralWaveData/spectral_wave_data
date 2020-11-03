# Build a simple example of C-program using
# the Fortran spectral-wave-data static library which is build first

cmake_minimum_required(VERSION 3.4.1)

set(PROG swd_c_app_example)
set(PROG_FFT swd_c_app_example_fft)
set(LIB spectral_wave_data)

project(${PROG} C Fortran)
project(${PROG_FFT} C Fortran)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_API_C ..)
set(DIR_SRC_C_APP ../../../api_examples)

include(${DIR_SRC_API_F}/Cmake/core_files.cmake)

# Core files plus kind values from the C directory
set(SRC_API_F
    ${SRC_CORE}
    ${DIR_SRC_API_C}/kind_values.f90 
    ${DIR_SRC_API_C}/spectral_wave_data.c 
)

set(SRC_C_APP
    ${DIR_SRC_C_APP}/application_swd_api.c
)

set(SRC_C_APP_FFT
    ${DIR_SRC_C_APP}/fft_example.c
)

add_library(${LIB} STATIC ${SRC_API_F})
target_include_directories(${LIB} INTERFACE ${DIR_SRC_API_C})

add_executable(${PROG} ${SRC_C_APP})
target_link_libraries(${PROG} PRIVATE ${LIB} fftw3)

add_executable(${PROG_FFT} ${SRC_C_APP_FFT})
target_link_libraries(${PROG_FFT} PRIVATE ${LIB} fftw3)