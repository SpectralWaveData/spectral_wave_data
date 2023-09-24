Changelog
=========

Version 1.0.0 (released 2023.09.XX)
"""""""""""""""""""""""""""""""""""

Minor New Features
------------------

-  The limit on file name length in constructor is removed. (Previous limit was 200 characters)
-  In Python, the constructor parameters `x0`, `y0`, `t0` and `beta` have default values `0.0`.
-  Support for Python :doc:`with-statement<api_specification_Python>` implemented for creating runtime context.
-  Default constructor added to the C++ implementation.
-  The python package provides a new function :doc:`airy.get_components()<tools>`
   for extracting the Airy spectral components from a `shp=6` swd file.
-  Support for swd files with only one time step.

Miscellaneous
-------------

-  "pip install spectral-wave-data" (from PyPI) is now supported on both Linux and Windows. 
   (Previous specific versions had to be compiled for actual Linux distro)
-  The Python *.whl (pip install) file on PyPI is generic for Python 3.x and Python 2.7. 
   (Previous PyPI version only accepted 2.7, 3.5, 3.6, 3.7 and 3.8) Python 2.7 will not be supported
   and will not work in the next release. 3.x releases not supported by the 
   `latest numpy <https://numpy.org/news/#releases>`_ may still work but will not be supported
   in future releases. Currently the latest numpy only supports python 3.9 and newer.
-  Update of demo C++ application code. Two examples are provided. One using C++ unique pointers.
   This is the reccommended method for applying C++ exception handling when calling the SWD constructor.
-  On Windows the distributed precompiled libraries are based on Visual Studio 2022 and 
   Intel ifort OneAPI 2023.1
-  General documentation is updated. 
-  Automatic API description from source code is improved.
-  More error detections implemented.

Version 1.0.0.rc1 (released 2020.04.02)
"""""""""""""""""""""""""""""""""""""""

First public release.

