# spectral_wave_data

An API for representation and interchange of ocean wave kinematics.

This Python package is part of the spectral_wave_data API developed by
the GitHub organization [SpectralWaveData](https://github.com/SpectralWaveData).

![Flow chart showing how SWD can be used](https://raw.githubusercontent.com/SpectralWaveData/spectral_wave_data/master/docs/source/figures/swd_scheme.png)


## Purpose and goals

The main goal of the SWD initiative is to provide an open API for boosting
research and industrial application of ocean wave kinematic.
Our goal is that by providing a simple, but powerful, interface we can tear
down walls between the fields of oceanography, marine hydrodynamics and
structural engineering. A motivating example: A specialist tool from the field
of oceanography can generate waves that can in turn be simulated in a fully
viscous Navier-Stokes calculation (CFD) to assess, e.g., the motion of a
floating wind turbine, without either of the tools knowing about each other.
Our hope is that SWD can become a universal language for representing waves.

The SWD API can be utilized by academic and commercial software developers 
to assess effects of ocean gravity waves due to the permissive MIT license
and the thorough documentation provided.


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
- Python 3 and 2.7

The spectral_wave_data Python package mainly contains a class **SpectralWaveData**
providing the official Python API. In addition, it provides a script **swd_meta**
for extracting typical meta data from a given SWD-file. Such files contain
the spectral information of the actual wave trains.

As explained in the documentation the class **SpectralWaveData** provides 
consistent exception handling.


## Usage

Typical usage in a Python script or program::

    from spectral_wave_data import SpectralWaveData
    
    with SpectralWaveData(file_swd="my.swd") as swd:

        swd.update_time(time=3.8143127)      # Maximum time can be obtained from: swd['tmax']
    
        zeta = swd.elev(x=5.3, y=12.4)
        acc_e = swd.acc_euler(x=7.3, y=-3.4, z=-8.3)
        acc_p = swd.acc_particle(x=7.3, y=-3.4, z=-8.3)
        phi2 = swd.grad_phi_2nd(x=7.3, y=-3.4, z=-8.3)
    
        print("Surface elevation at (x, y) = ", zeta)
        print("Euler (local) acceleration in x-direction at (x,y,z) = ", acc_e.x)
        print("Particle acceleration in x-direction at (x,y,z) = ", acc_p.x)
        print("d^2(potential)/dxdz at (x,y,z) = ", phi2.xz)

Meta data from an actual SWD file can be extracted using the console script swd_meta::

    > swd_meta my.swd
    version: 1.0.0
    prog:    raschii-1.0.4
    date:    2020:10:22 18:37:55
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

To install spectral_wave_data from PyPI execute::

    pip install spectral_wave_data
    
in your actual Python environment on Linux or Windows. 
Python 2.7 and 3.x are supported.

On Linux spectral_wave_data is compiled on the manylinux2014 platform to confirm
with most modern Linux systems.

For other systems you need to build the package from source.
To compile this Python package from source you need access to a 
recent Fortran compiler. E.g gfortran-10 or Intel Fortran 2019 or later.

Due to the complexity of building from source we refer to the 
[documentation](https://spectral-wave-data.readthedocs.io/) and
the [spectral_wave_data](https://github.com/SpectralWaveData/spectral_wave_data)
GitHub repository for related instructions and download of the source code.


## Copyright and license

This software is copyrighted and licensed open-source under the MIT License.

See the [spectral_wave_data](https://github.com/SpectralWaveData/spectral_wave_data)
GitHub repository for details.
