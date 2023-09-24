# Build SpectralWaveData library for Python package

cmake_minimum_required(VERSION 3.20)

set(LIB SpectralWaveData)

project(${LIB} C Fortran)

if(WIN32)
    set_property (DIRECTORY PROPERTY VS_STARTUP_PROJECT ${LIB})
    # To static link Intel libraries into SWD (not yet recognized by VS, to be set manually in VS. see Readme file)
    set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS_RELEASE} /libs:static")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_BUILD_TYPE Release)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_API_C ../../c)
set(DIR_SRC_API_PY ..)

include(${DIR_SRC_API_F}/Cmake/core_files.cmake)

# Core files plus kind values from the Python directory
set(SRC_API_F
    ${SRC_CORE}
    ${DIR_SRC_API_PY}/kind_values.f90
)

add_library(${LIB} SHARED ${SRC_API_F})
target_include_directories(${LIB} INTERFACE ${DIR_SRC_API_C})

