*****
Tools
*****

In this section we describe some auxiliary tools related to **spectral_wave_data**.

=======
swd2vtk
=======

.. figure:: figures/swd2vtk_64x64.png
   :align: left

A program **swd2vtk** is available for advanced three-dimensional visualization of SWD wave fields.
For a user specified domain a set of VTK-files is produced. The user selects
a set of scalar and vector fields to be generated. VTK-files can be animated or further post-processed
using the open-source program `ParaView <https://www.paraview.org/>`_ or other tools. More details are
described in the :doc:`swd2vtk user's guide<tools_swd2vtk>`

Binary distributions of **swd2vtk**, for Windows and Linux, can be downloaded from the release tab of the
GitHub repository **spectral_wave_data**. The binaries are distributed under the MIT-license.
Hence, you can basically apply **swd2vtk** as you like for free.

.. toctree::
   :hidden:

   tools_swd2vtk

==========
Airy Waves
==========

A general set of Airy waves can be written directly to a SWD file (:doc:`shape 6<shape_6>`) using
the :mod:`tools.airy` module included in **spectral_wave_data**.

In this example we apply two wave components.

.. code-block:: python

    from spectral_wave_data.tools import airy

    airy.write_swd('my_airy.swd',
                   amps=[2.0, 3.0], dirs=[180.0, 70.0], phases=[50.0, 0.0], Twaves=[10.5, 13.0],
                   depth=32.0, grav=9.81, is_deg_dirs=True, is_deg_phases=True)

The :doc:`documentation <tools_airy_write_swd>` of this function describes relevant
optional parameters and further examples.

The spectral components stored in an SWD file (shp=6) can be extracted into a dictionary 
using the syntax.

.. code-block:: python

    from spectral_wave_data.tools import airy

    >>> res = airy.get_components('my_airy.swd')
    >>> res
    {"n" : 2, 
     "depth" : 32.0,
     "grav" : 9.81,
     "amps" : [2.0, 3.0], 
     "dirs_deg" : [180.0, 70.0], 
     "phases_deg" : [50.0, 0.0], 
     "Twaves" : [10.5, 13.0],
     "kwaves" : [0.041875288176738856, 0.031260688477796145],
     "omegas" : [0.5983986006837702, 0.483321946706122],
     "Lwaves" : [150.04518370502365, 200.99318387189103]]}

.. toctree::
   :hidden:

   tools_airy_write_swd

===================
Fenton-Stream Waves
===================

The `raschii <https://pypi.org/project/raschii/>`_ Python package for making non-linear
regular waves supports output of SWD files. An example of producing a Fenton-Stream wave
is given by this code snippet.

.. code-block:: python

    import raschii

    WaveModel, AirModel = raschii.get_wave_model('Fenton')
    wave = WaveModel(height=5.0, depth=15.0, length=200.0, N=30)
    wave.write_swd('my_fenton.swd', dt=0.05, nperiods=50)

This method applies the :mod:`shp=2` shape class. The module :file:`swd_tools.py` in the raschii
distribution demonstrates how to create SWD files directly from Python.

Raschii also has an online graphical `wave calculator <https://raschii.readthedocs.io/en/latest/raschii_dart.html>`_.

============
Stokes Waves
============

The `raschii <https://pypi.org/project/raschii/>`_ Python package may also output Stokes
waves in the SWD format.

.. code-block:: python

    import raschii

    WaveModel, AirModel = raschii.get_wave_model('Stokes')
    wave = WaveModel(height=5.0, depth=15.0, length=200.0, N=5)
    wave.write_swd('my_stokes.swd', dt=0.05, nperiods=50)

This method applies the :mod:`shp=2` shape class.

=====
WAMOD
=====

.. figure:: figures/Wamod.png
   :align: left

DNV has
`released the program WAMOD <https://www.dnv.com/services/hydrodynamic-analysis-of-non-linear-waves-in-random-sea-states-wamod-177704>`_
for simulation of irregular seas based on HOS(M), the Higher-Order-Spectral-Method.
The typical input is a wave spectrum and the output is a SWD file providing linear or nonlinear
wave fields. Long and short crested seas can be simulated in finite or infinite water depth.

WASIM can apply the SWD API for incident wave kinematics and Froude-Krylov forces.
There are no restrictions regarding which wave generator produced the SWD file.
WASIM is a Rankine panel program for simulating wave induced responses of marine structures.

WAMOD and WASIM have commercial licenses and support.

===============
Other programs
===============

Other software utilizing SWD will be added to this page when they are released.
Please contact the maintainers when new software appear.



