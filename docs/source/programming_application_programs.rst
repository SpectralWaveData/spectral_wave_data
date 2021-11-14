.. _ref_programming_applications:

====================
Application programs
====================

There are no restrictions on programming styles applied in
:ref:`application programs<application-program>`. However, the SWD API is object oriented
in order to simplify its application.
An instance of a specific class takes care of
all related kinematics and swd-file operations.
The application programmer do not consider I/O operations at all.

The API interface to be applied in :ref:`application programs<application-program>`
is basically the very same for all of the following programming languages:
Fortran-2008, C/C++ and Python-2/3. Fully working :ref:`application programs<application-program>`
demonstrating most features of the API in all of these languages follows

.. toctree::

   application_swd_app_python
   application_swd_app_Cpp
   application_swd_app_C
   application_swd_app_fortran

---------------------
Compiling and linking
---------------------

For using precompiled binaries checkout the PREREQUISITES.md file in the root of
the GitHub repository.

Due to applications of object oriented features, only modern Fortran compilers are able
to compile the source code. Any ISO standard conforming Fortran-2008 compiler can be applied.
In practice you need a very recent Fortran compiler to successfully compile this package. It will
successfully compile using gfortran 9.x and Intel 2019/2020 compilers. Older compilers will
not work due to lack of ISO standard support.

Regarding C++, C and Fortran applications the source code from this repository is written
in ISO-standard Fortran-2008, C and C++. Consequently you may compile and link the source from this
repository directly with your application if you apply standard conforming compilers.
There are `Cmake <https://cmake.org>`_ files to demonstrate how to do this for the test
applications listed above. These scripts have been tested on Windows-10 and Linux-Ubuntu.

.. note::
  When compiling for C/C++ these compilers must be binary compatible with the actual
  Fortran compiler. E.g. If you apply GNU compilers for C/C++ you should apply gfortran.
  If you apply Microsoft/Intel C/C++ compilers you should apply Intel Fortran.

.. note::
  For gfortran you need gfortran-10
  if you need the possibility to have several concurrent SWD objects connected to the same SWD file.
  For Intel compilers this is not an issue.
