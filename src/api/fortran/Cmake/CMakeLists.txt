# Build spectral_wave_data Fortran library
# The build a small test program using the library

cmake_minimum_required(VERSION 3.4.1)

set(PROG swd_f_app_example)
set(LIB spectral_wave_data)

project(${PROG} Fortran)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_APP ../../../api_examples)

set(SRC_SWDLIB
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
  ${DIR_SRC_API_F}/swd_write_shape_1_or_2.f90
  ${DIR_SRC_API_F}/swd_write_shape_3.f90
  ${DIR_SRC_API_F}/swd_write_shape_4_or_5.f90
  ${DIR_SRC_API_F}/swd_write_shape_6.f90
  ${DIR_SRC_API_F}/kind_values.f90
  ${DIR_SRC_API_F}/swd_version.f90)

set(SRC_APP
    ${DIR_SRC_APP}/application_swd_api.f90
)

add_library(${LIB} ${LIBTYPE} ${SRC_SWDLIB})
add_executable(${PROG} ${SRC_APP})
target_link_libraries(${PROG} PRIVATE ${LIB}) 
