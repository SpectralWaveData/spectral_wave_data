# All the Fortran implementation source files except kind_values.f90
# You need to define DIR_SRC_API_F before include()-ing this file
set(SRC_CORE
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
  ${DIR_SRC_API_F}/swd_write_shape_1_or_2.f90
  ${DIR_SRC_API_F}/swd_write_shape_3.f90
  ${DIR_SRC_API_F}/swd_write_shape_4_or_5.f90
  ${DIR_SRC_API_F}/swd_write_shape_6.f90
  ${DIR_SRC_API_F}/swd_version.f90)

# Bundle the Intel compiler libraries when compiling shared libraries
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (UNIX)
        message("SWD: enabling -static-intel for .so files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -static-intel")
    endif()
    if (WIN32)
        message("SWD: enabling /static for .lib files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} /libs:static /threads")
    endif()
endif()


# Bundle the Intel compiler libraries when compiling shared libraries on Linux
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND UNIX)
#    message("SWD: enabling -static-intel for .so files")
#    set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -static-intel")
#endif()