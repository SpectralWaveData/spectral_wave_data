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
  ${DIR_SRC_API_F}/swd_version.f90
  ${DIR_SRC_API_F}/swd_fft.f90
  ${DIR_SRC_API_F}/swd_fft_fftw3.f90)

# Bundle the Intel compiler libraries when compiling shared libraries
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    get_filename_component(COMP_PATH ${CMAKE_Fortran_COMPILER} DIRECTORY)
    set(MKLROOT "${COMP_PATH}/../mkl")
    message(STATUS "MKLROOT: ${MKLROOT}")
    INCLUDE_DIRECTORIES(${MKLROOT}/include/fftw/)
    message(STATUS "FFTW_INCLUDE_DIRS: ${MKLROOT}/include/fftw/")
    if (UNIX)
        message("SWD: enabling -static-intel for .so files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -static-intel")
    endif()
    if (WIN32)
        message("SWD: enabling /static for .lib files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib /libs:static /threads")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND UNIX)
    find_path (FFTW_INCLUDE_DIRS fftw3.f03)
    message(STATUS "FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
endif()


# Bundle the Intel compiler libraries when compiling shared libraries on Linux
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND UNIX)
#    message("SWD: enabling -static-intel for .so files")
#    set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -static-intel")
#endif()
