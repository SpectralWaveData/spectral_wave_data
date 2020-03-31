************************
Verification and Testing
************************

It is important to know that for **spectral_wave_data** the concept of verification
is much more relevant than validation.

The quality of a specific wave field should only depend on the :ref:`wave generator<wave-generator>`.
Hence, **spectral_wave_data** should only reproduce the kinematics as it is defined by the
contents of the SWD file and the :doc:`mathematical spectral definition <theory>`.
However, this reproduction should be of high accuracy.

--------------------------
The verification procedure
--------------------------

Each of the temporal interpolation schemes described in the :doc:`theory <theory>`
section is consistent to a certain order.
Consequently, for all spectral amplitudes defined by polynomials of the actual order,
these schemes should produce exact temporal evaluations.

For each temporal interpolation scheme, random polynomial spectral amplitudes of the corresponding order
are constructed and stored in SWD files.
Consequently, the associated waves have no physical relevance.
Wave kinematics are then evaluated using two different procedures:

1.  The Python **spectral_wave_data** package, which apply the Fortran and C implementations,
    evaluates all supported kinematics from the constructed SWD files.
2.  A python test script combines the analytical polynomial spectral amplitudes, with the actual mathematical
    spectral definition. By using symbolic mathematics
    (`SymPy <https://www.sympy.org>`_) all kinematics are evaluated analytically using symbolic differentiation of
    the formulas for wave potential and surface elevation. Consequently, the explicit formulas in the
    :doc:`theory <theory>` sections for e.g. particle velocities are not applied in these evaluations.
    The symbolic calculations are all evaluated directly in the application specific time and space
    coordinate system.

Using `pytest <https://docs.pytest.org>`_ all kinematics from both the **spectral_wave_data** package and
the symbolic computations described above are calculated and compared. Any deviations are reported.

This verification procedures is applied for all spectral shape definitions and supported implementation
schemes, all temporal interpolation schemes, and a representative set of different constructor parameters
(x0, y0, t0 and beta).

For shape 6 (Airy waves) the auxilary tool included in **spectral_wave_data** is
applied when writing the relevant test SWD files.
Hence, :meth:`spectral_wave_data.tools.airy.write_swd` is verified too.

The master/origin branch of this repository does not update unless all these verification tests
pass successfully.

-------
Testing
-------

All tests are executed using `pytest <https://docs.pytest.org>`_ and scripts located
in the :file:`tests/python` directory of the **spectral_wave_data** repository.

In this directory there are some relevant scripts::

  install_requirements_tests.bat
  run_all_tests.bat
  run_failed_test.bat

There are similar scripts for Linux. The script :file:`install_requirements_tests.bat`  will if needed,
install all required Python packages (from `PyPI <https://pypi.org/>`_) for running the tests and
creating the test report.

Running :file:`run_all_tests.bat` will run all tests and display a html report when finished. Normal output
from the package and scripts are suppressed by `pytest <https://docs.pytest.org>`_.

The script :file:`run_failed_test.bat` reruns only the first failed test and provide detailed output including
local variables from the python scripts.

Running all tests may take a long time but is required in order to submit
a pull request for updating this GitHub repository.

New features, like new spectral formulations, will not be included in the master/origin branch before proper
`pytest <https://docs.pytest.org>`_-scripts have been established. If needed, we may assist...


.. note::

  The test suite described in this section apply test features not supported in Python 2.7. Hence it
  is only possible to run these tests in Python 3.5 and newer.

  The core package **spectral_wave_data** seems to run well on Python-2.7 too, but support for 2.x may be
  dropped in the future if dependency packages drop support for 2.x.
  This is a general trend in
  the Python community because January 1st, 2020, Python 2.7 officially reached the end of
  life and will no longer receive security updates, bug fixes, or other improvements going forward.

  As a consequence we strongly recommend to apply Python-3 for future proof utilization
  of  **spectral_wave_data**.