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

include(${DIR_SRC_API_F}/Cmake/core_files.cmake)

# Core files plus kind values from the C++ directory
set(SRC_API_F
    ${SRC_CORE}
    ${DIR_SRC_API_CPP}/kind_values.f90
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
