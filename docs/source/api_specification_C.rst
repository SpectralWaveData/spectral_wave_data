**
C
**

C does not support object-oriented-programming. However, an object based
approach is applied. All functions apply the prefix naming convention :attr:`swd_api_*`.
The API is defined in the header file :doc:`spectral_wave_data.h<c_spectral_wave_data_h>`.

The data types
:attr:`real_swd`, :attr:`vector_swd`, :attr:`vector_2nd_phi_swd` and :attr:`vector_2nd_elev_swd`
are also defined in the C-header file.

If your application apply float (not double) for the SWD-API you need to:

1.  Set the macro :attr:`SWD_API_FLOAT` before the header file is included.
2.  The fortran library should be compiled with `kind_swd_c = c_float` in
    the file :file:`kind_values.f90`.

By default it is assumed that the C-interface apply `double`. Hence `kind_swd_c = c_double` in
the file :file:`kind_values.f90` and the macro :attr:`SWD_API_FLOAT` should be unset.

.. toctree::
   :hidden:

   c_spectral_wave_data_h.rst

--------------
Implementation
--------------

In the current version the implementation of these functions are obtained
from the fortran implementation using the standard ISO_C_BINDINGS defined in
:doc:`spectral_wave_data_c.f90<fortran_swd_c>`.

.. toctree::
   :hidden:

   fortran_swd_c.rst
