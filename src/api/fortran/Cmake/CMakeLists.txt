# Build spectral_wave_data Fortran library
# The build a small test program using the library

cmake_minimum_required(VERSION 3.4.1)

set(PROG swd_f_app_example)
set(LIB spectral_wave_data)

project(${PROG} Fortran)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_APP ../../../api_examples)

include(core_files.cmake)

# Core files plus kind_values.f90
set(SRC_SWDLIB
  ${SRC_CORE}
  ${DIR_SRC_API_F}/kind_values.f90)

set(SRC_APP
    ${DIR_SRC_APP}/application_swd_api.f90
)

add_library(${LIB} ${LIBTYPE} ${SRC_SWDLIB})
add_executable(${PROG} ${SRC_APP})
target_link_libraries(${PROG} PRIVATE ${LIB}) 
