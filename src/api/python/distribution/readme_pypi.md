# spectral_wave_data

An API for ocean wave kinematics

This Python package is part of the spectral_wave_data API developed by
the GitHub organization [SpectralWaveData](https://github.com/SpectralWaveData).

![alt text](https://github.com/SpectralWaveData/spectral_wave_data/blob/master/docs/source/figures/repository_readme_600x487.png)


## Purpose and goals

The main idea is to provide an open API for boosting research and 
industrial application of spectral ocean wave kinematic.
This API can be utilized by academic and commercial software developers 
to assess effects of ocean gravity waves. Not limited to spectral kinematics, 
but also, as advanced boundary conditions for more sophisticated wave 
models and to study related effects on marine structures.

By using this simple but powerful interface another goal with this API 
is to tear down walls between the fields of oceanography, marine 
hydrodynamics (flow around marine structures) and structural engineering.


## What kind of waves?

There are many spectral formulations out there. Our goal is not to 
implement all of them, but to make the API generic in the sense that 
most formulations can apply this API, and that
additional formulations can easily be tested or added at a later stage
without breaking the current API. At this stage we have implemented 
long and short crested waves propagating in infinite or constant water 
depth. For long crested seas support for varying depth is also included.

The current version has been used for describing, Airy, Stokes, 
Stream waves and Higher-Order-Spectral-Method waves of arbitrary order.
Detailed kinematics of simulated rogue waves can be evaluated.


## Documentation

Detailed [documentation](https://spectral-wave-data.readthedocs.io/)
describing theory, examples, implementation and related 
tools is hosted on [ReadTheDocs](https://readthedocs.org/).


## The software

This is the python package. However, in general 
[spectral_wave_data](https://github.com/SpectralWaveData/spectral_wave_data)
provides API for several other programming languages:

- ISO Fortran-2008
- ISO C / C++
- Python 2 and 3

The spectral_wave_data Python package mainly contains a class **SpectralWaveData**
providing the official Python API. In addition, it provides a script **swd_meta**
for extracting typical meta data from a given SWD-file. Such files contain
the spectral information of the actual wave trains.

As explained in the documentation the class **SpectralWaveData** provides 
consistent exception handling.


## Usage

Typical usage in a Python script or program::

    from spectral_wave_data import SpectralWaveData
    
    # x0/y0/t0/beta are parameters shifting the spatial and temporal reference frame
    swd = SpectralWaveData(file_swd="my.swd", x0=0.0, y0=0.0, t0=0.0, beta=180.0)

    swd.update_time(time=3.8143127)      # Maximum time can be obtained from: swd['tmax']
    
    zeta = swd.elev(x=5.3, y=12.4)
    acc_e = swd.acc_euler(x=7.3, y=-3.4, z=-8.3)
    acc_p = swd.acc_particle(x=7.3, y=-3.4, z=-8.3)
    phi2 = grad_phi_2nd(x=7.3, y=-3.4, z=-8.3)
    
    print("Surface elevation at (x, y) = ", zeta)
    print("Euler (local) acceleration in x-direction at (x,y,z) = ", acc_e.x)
    print("Particle acceleration in x-direction at (x,y,z) = ", acc_p.x)
    print("d^2(potential)/dxdy at (x,y,z) = ", phi2.xy)

Meta data from an actual SWD file can be extracted using the console script swd_meta::

    >>> swd_meta my.swd
    version: 1.0.0
    prog:    raschii-1.0.3
    date:    2020:01:22 19:57:55
    fmt:     100
    shp:     2
    amp:     1
    tmax:    6.3000000938773155
    dt:      0.10000000149011612
    nsteps:  64
    nstrip:  0
    order:   -1
    d:       32.0
    n:       50
    sizex:   220.00000561733003
    lmax:    220.00000561733003
    lmin:    4.400000112346601
    dk:      0.028559932485222816
    cid:     {'model': 'Fenton', 'T': 12.792885811907514, 'height': 18.5, ...}

Output and further examples are explained in the 
[online documentation](https://spectral-wave-data.readthedocs.io/).


## Installation

To install spectral_wave_data is usually obtained executing::

    pip install spectral_wave_data
    
in your actual Python environment on Linux or Windows-10. 
Python 2.7 or 3.5, 3.6, 3.7 and 3.8 are supported.

On windows you may have to download and install the
[Intel redistributable library for Fortran 2020](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2020-compilers-for-windows)
using default options, unless you have a recent Intel Fortran compiler installed.

On Linux spectral_wave_data is compiled using gfortran. The related GNU libraries 
should most likely already be available in your system.

For other systems you need to build the package from source.
To compile this Python package from source you need access to a 
recent Fortran compiler. E.g gfortran-10 or Intel Fortran 2019 or later.

Due to the complexity of building from source we refer to the 
[documentation](https://spectral-wave-data.readthedocs.io/) and
the [spectral_wave_data](https://github.com/SpectralWaveData/spectral_wave_data)
GitHub repository for related instructions and download of the source code.

### Installation problems?

In some rare and odd Python environments users have encountered issues loading
the spectral_wave_data in combinations with e.g. [numpy](https://pypi.org/project/numpy/).
The best solution is to reinstall your Python system in a proper way.
However, a quick-fix is usually to load spectral_wave_data before numpy because
spectral_wave_data requires a more recent version of some Intel DLLs on Windows.

This rare issue is not present if you apply recommended procedures for 
setting up your Python system. Hence, consistently using e.g. pyenv or Anaconda environments.
In general, manual fumbling of system paths is always a source of subtle problems
and should be avoided.


## Copyright and license

This software is copyrighted and licensed open-source under the MIT License.

See the [spectral_wave_data](https://github.com/SpectralWaveData/spectral_wave_data)
GitHub repository for details.
