***
C++
***

The official C++ API is defined in the header file
:doc:`SpectralWaveData.h<cpp_SpectralWaveData_h>`.
It defines the generic class :class:`SpectralWaveData` to be applied in applications.
All member functions and input arguments are explained in that header file.

The data types
:attr:`real_swd`, :attr:`vector_swd`, :attr:`vector_2nd_phi_swd` and :attr:`vector_2nd_elev_swd`
are defined in the C-header file :doc:`spectral_wave_data.h<c_spectral_wave_data_h>`.

If your C++ application apply float (not double) for the SWD API interface you need to:

1.  Set the macro :attr:`SWD_API_FLOAT` before the header file is included.
2.  The fortran library should be compiled with `kind_swd_c = c_float` in
    the file :file:`kind_values.f90`.

By default it is assumed that the C/C++ interface apply `double`. Hence `kind_swd_c = c_double` in
the file :file:`kind_values.f90` and the macro :attr:`SWD_API_FLOAT` should be unset.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - C++ Exceptions
     - Usage
   * - :exc:`SwdException`
     - Base class for all SWD specific exceptions. It inherits from :exc:`std::runtime_error`.
   * - :exc:`SwdFileCantOpenException`
     - Typical if SWD file is not an existing file.
   * - :exc:`SwdFileBinaryException`
     - SWD file does not apply float/little-endian
   * - :exc:`SwdFileDataException`
     - Error during reading and checking SWD data
   * - :exc:`SwdInputValueException`
     - Input arguments for class methods are not sound
   * - :exc:`SwdAllocateException`
     - Not able to allocate internal SWD storage

.. toctree::
   :hidden:

   cpp_SpectralWaveData_h.rst
   cpp_SpectralWaveData_cpp.rst

--------------
Implementation
--------------

In the current version the implementation of the class :class:`SpectralWaveData`
is obtained by wrapping the C-implementation as done in the file
:doc:`SpectralWaveData.cpp<cpp_SpectralWaveData_cpp>`.

