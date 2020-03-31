***************************************
spectral_wave_data.tools.airy.write_swd
***************************************

The following function available in the sub-package :mod:`spectral_wave_data.tools`
generates a SWD-file for a general set of Airy waves. All conventions follow the
definitions described by :doc:`shape 6<shape_6>`.

.. autofunction:: spectral_wave_data.tools.airy.write_swd

--------
Examples
--------

First we establish access to the function:

.. code-block:: python

    from spectral_wave_data.tools import airy

In the first example we assume that wave numbers are explicitly defined and the water depth is infinite.
In all examples we assume that both phases and propagation directions are given in the unit [deg]

.. code-block:: python

    amps = [...]   # Single wave amplitudes for each component
    dirs = [...]   # Wave propagation directions for each component
    phases = [...] # Phase angles for each component
    kwaves = [...] # Wave numbers for each component
    airy.write_swd('my_airy.swd', amps, dirs, phases, kwaves=kwaves,
                   depth=-1.0, is_deg_dirs=True, is_deg_phases=True)

In case wave numbers are implicitly defined by wave periods:

.. code-block:: python

    Twaves = [...] # Wave periods for each component
    airy.write_swd('my_airy.swd', amps, dirs, phases, Twaves=Twaves,
                   depth=-1.0, is_deg_dirs=True, is_deg_phases=True)

In case wave numbers are implicitly defined by wave frequencies:

.. code-block:: python

    omegas = [...] # Wave frequencies for each component
    airy.write_swd('my_airy.swd', amps, dirs, phases, omegas=omegas,
                   depth=-1.0, is_deg_dirs=True, is_deg_phases=True)

In case wave numbers are implicitly defined by wave lengths:

.. code-block:: python

    Lwaves = [...] # Wave lengths for each component
    airy.write_swd('my_airy.swd', amps, dirs, phases, Lwaves=Lwaves,
                   depth=-1.0, is_deg_dirs=True, is_deg_phases=True)

In case `Twaves` or `omegas` is applied the linear wave dispersion relation is applied to establish
the corresponding wave numbers to be stored in the SWD file. This function is available from the
same Python module:

.. autofunction:: spectral_wave_data.tools.airy.omega2kwave

