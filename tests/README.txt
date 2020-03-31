This directory contains all test codes related to spectral-wave-data

The tests are sorted into language dependent subdirectories:

python
------

The main testing of spectral-wave-data is carried out in python using pytest.
If these tests works than you can also conclude that the fortran and C code is working too
because the Python implementation wraps the C implementation that wraps the Fortran implementation.

fortran
-------

Due to the extensive python testing, only a few explicit fortran tests are present.
These tests are manual executed based on individual Cmake configurations.

C
-

See testing of the python interface which apply the C-interface.


C++
---

No explicit tests yet. The C++ implementations is just a wrapper of the C-implementation.
If the spectral-wave-data\src\api_examples\application_swd_api.cpp works using 
spectral-wave-data\src\api\cpp\Cmake the complete interface works.


