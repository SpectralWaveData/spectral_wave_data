**************
Python 2 and 3
**************

The Python package **spectral_wave_data** contains a module
:doc:`spectral_wave_data.py<python_spectral_wave_data>`.
It defines the generic class :class:`SpectralWaveData`
to be applied in applications based on Python-2 or Python-3.

-----------------------
Constructor and methods
-----------------------

Detailed documentation of class members is provided :doc:`here <api_specification_Python_methods>`.

.. toctree::
   :hidden:
   :caption: autdoc_python

   api_specification_Python_methods

------------------
Exception handling
------------------

Related to the class :class:`SpectralWaveData` is a set of class specific exceptions
defined in the same module. These exceptions as described in the table below, can be
imported and applied in application programs and scripts.

>>> from spectral_wave_data import SpectralWaveData, SwdError, SwdInputValueError, ...
>>> ...
>>> try:
>>>    swd = SpectralWaveData('my_waves.swd', x0=0.0, y0=0.0, t0=0.0, beta=180.0)
>>> except SwdError as err:
>>>    print(err)  # Print actual error message
>>>    # You may do some application specific error recovery before the show must go on...
>>> ...

If no try/except block is applied, the application program will abort and
print the exception and backtrace as a normal Python crash.

The only other associated class methods that may throw exceptions are :meth:`update_time`,
:meth:`convergence`, :meth:`strip` and :meth:`get`.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Python Exceptions
     - Usage
   * - :exc:`SwdError`
     - Base class for all SWD specific exceptions
   * - :exc:`SwdFileCantOpenError`
     - Typical if SWD file is not an existing file.
   * - :exc:`SwdFileBinaryError`
     - SWD file does not apply float/little-endian
   * - :exc:`SwdFileDataError`
     - Error during reading and checking SWD data
   * - :exc:`SwdInputValueError`
     - Input arguments for class methods are not sound
   * - :exc:`SwdAllocateError`
     - Not able to allocate internal SWD storage

--------
Metadata
--------

To extract metadata from a SWD file you may apply the method :mod:`swd.get(key)`,
or using the more Pythonesque syntax :mod:`swd[key]`, where `key` is a string to identify
the requested metadata. A key is either:

  * A relevant parameter from the :doc:`SWD-file format description <swd_format>`.
  * A constructor parameter.
  * A key from the table below.

Only scalar meta data is supported in this version.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Key
     - Returned value
   * - version
     - | The repository version number of **spectral_wave_data** applied in this
       | Python distribution.
   * - class
     - Name of the specialized SWD-class handling this object. (Fortran class)
   * - tmax
     - | Maximum allowed application time consistent with the header of
       | the SWD file and the constructor parameter `t0=0.0`. [s]
   * - lmin
     - Shortest wave length component. [m]
   * - lmax
     - Longest wave length component. [m]
   * - sizex
     - Periodic length of wave domain in x-direction (swd-system). [m]
   * - sizey
     - Periodic length of wave domain in y-direction (swd-system). [m]

A  :exc:`SwdInputValueError` exception is raised if the actual key is not relevant
for the actual SWD class.

^^^^^^^^^^^^^^^^^^^^^^^^^^^
The script swd_meta
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For convenience this Python wheel package includes the script :file:`swd_meta`
listing the relevant metadata for a given SWD-file. It runs on Windows and Linux.

In a terminal window with access to your installed **spectral_wave_data** package you can
invoke it like in this example where :file:`my.swd` is the name of the actual SWD-file:

>>> swd_meta my.swd
version: 1.0.0-beta.9
prog:    raschii-1.0.3.dev0
date:    2020:01:22 19:57:55
fmt:     100
shp:     2
amp:     1
tmax:    6.3000000938773155
dt:      0.10000000149011612
nsteps:  64
nstrip:  0
order:   -1
depth:   32.0
n:       50
sizex:   220.00000561733003
lmax:    220.00000561733003
lmin:    4.400000112346601
dk:      0.028559932485222816
cid:     {'model': 'Fenton', 'T': 12.792885811907514, 'height': 18.5, 'depth': 32.0, 'N': 50, 'air': 'NoneType', 'g': 9.81, 'c': 17.19705805512825, 'relax': 0.5}


--------------
Implementation
--------------

In the current version the implementation of the class :class:`SpectralWaveData`
is obtained by wrapping the C-implementation using the standard Python module
:mod:`ctypes` in the package module :doc:`swd_c_interface.py<python_swd_c_interface>`.

.. toctree::
   :hidden:

   python_spectral_wave_data.rst
   python_swd_c_interface.rst
