# -*- coding: utf-8 -*-
    
"""
:platform: Linux, Windows, python 2.7 and 3.x
:synopsis: Create an SWD file using shape = 6 (Airy waves)

Author  - Jens Bloch Helmers, DNVGL
Created - 2019-11-17
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

from struct import pack, unpack
from datetime import datetime
import numpy as np
import sys

__all__ = ['write_swd', 'omega2kwave', 'get_components']


def write_swd(path, amps, dirs, phases, kwaves=None, Twaves=None,
              omegas=None, Lwaves=None, depth=-1.0, grav=9.81,
              Lscale=1.0, is_deg_dirs=True, is_deg_phases=True):
    """Creates a SWD file of Airy wave components (shape 6).

    One and only one of the sequences `kwaves`, `Twaves`,
    `omegas` or `Lwaves` should be specified.

    Parameters
    ----------
    path : str
        The path of the new SWD file
    amps : sequence of floats
        Sequence of wave amplitudes. [m] (single amplitudes)
    dirs: sequence of floats
        Sequence of wave propagation directions.
        The input unit is [deg] or [rad] depending on the value of `is_deg_dirs`.
    phases: sequence of floats
        List of wave phases as defined in shape 6.
        The input unit is [deg] or [rad] depending on the value of `is_deg_phases`.
    Twaves: Sequence of floats, optional
        Sequence of wave periods [sec]
    kwaves: Sequence of floats, optional
        Sequence of wave numbers [1/m]
    omegas: Sequence of floats, optional
        Sequence of wave frequencies [rad/sec]
    Lwaves: Sequence of floats, optional
        Sequence of wave lengths [m]
    depth: float, optional
        Actual water depth [m]  depth < 0 indicates infinite water depth. (default)
    grav: float, optional
        Acceleration of gravity in wave generator. [m/s^2] (ref item in SWD-file)
    Lscale: float, optional
        Number of length units per meter in wave generator. [-] (ref item in SWD-file)
    is_deg_dirs: bool, optional
        True if dirs are given in the unit deg, (rad if False)
    is_deg_phases: bool, optional
        True if phases are given in the unit deg, (rad if False)

    Returns
    -------
    None
       Nothing

    Raises
    ------
    IOError
        Path can not be opened for writing the SWD file.
    AssertionError
        Inconsistent input. E.g. length of sequences differ.

    """

    prog = 'spectral_wave_data.tools.airy'
    dtim = datetime.now().strftime("%Y:%m:%d %H:%M:%S")
    wave_generator_data = str({'path': path,
                               'amps': amps,
                               'dirs': dirs,
                               'phases': phases,
                               'kwaves': kwaves,
                               'omegas': omegas,
                               'Twaves': Twaves,
                               'Lwaves': Lwaves,
                               'depth' : depth,
                               'grav': grav,
                               'Lscale': Lscale,
                               'is_deg_dirs': is_deg_dirs,
                               'is_deg_phases': is_deg_phases,
                               })

    nwaves = len(amps)
    # Convert to wave numbers if not given...
    if kwaves is None:
        if Twaves is not None:
            omegas = 2.0 * np.pi / np.array(Twaves)
        if omegas is None:
            kwaves = 2.0 * np.pi / np.array(Lwaves)
        else:
            kwaves = np.array([omega2kwave(w, depth, grav) for w in omegas])

    assert len(dirs) == len(phases) == len(kwaves) == nwaves

    if is_deg_dirs:
        dirs = np.array(dirs) * (np.pi / 180.0)
    if is_deg_phases:
        phases = np.array(phases) * (np.pi / 180.0)

    out = open(path, 'wb')
    out.write(pack('<f', 37.0221))  # Magic number
    out.write(pack('<i', 100))  # fmt
    out.write(pack('<i', 6))  # shp for long-crested constant finite depth
    out.write(pack('<i', 1))  # amp (both elevation and field data)
    out.write(pack('<30s', prog.ljust(30).encode('utf-8')))  # prog name
    out.write(pack('<20s', dtim.ljust(20).encode('utf-8')))  # date
    nid = len(wave_generator_data)
    out.write(pack('<i', nid))  # length of input file
    out.write(pack('<{0}s'.format(nid), wave_generator_data.encode('utf-8')))  # Input file
    out.write(pack('<f', grav))  # acc. of gravity
    out.write(pack('<f', Lscale))  # lscale
    out.write(pack('<i', 0))  # nstrip
    out.write(pack('<i', 0))  # nsteps
    out.write(pack('<f', -1.0)) # dt
    out.write(pack('<i', 0))  # order
    out.write(pack('<i', nwaves))
    out.write(pack('<f', depth))
    for i in range(nwaves):
        out.write(pack('<f', amps[i]))
        out.write(pack('<f', kwaves[i]))
        out.write(pack('<f', dirs[i]))
        out.write(pack('<f', phases[i]))
    out.close()


def get_components(file_swd):
    """Return the Airy spectral components stored in file_swd

    Parameters
    ----------
    file_swd : str
        The path of the exisiting SWD file

    Returns
    -------
    res : dictionary with the following keys:

    'n' : Number of spectral components
    'depth' : Actual water depth [m]  depth < 0 indicates infinite water depth.
    'grav' : Acceleration of gravity [m/s^2]
    'amps' : Sequence of wave amplitudes. [m] (single amplitudes)
    'dirs_deg' : Sequence of wave propagation directions as defined in shape 6. (deg)
    'phases_deg' : List of wave phases as defined in shape 6. (deg)
    'Twaves' : Sequence of wave periods [sec]
    'Lwaves' : Sequence of wave lengths [m]
    'kwaves' : Sequence of wave numbers [1/m]
    'omegas' : Sequence of wave frequencies [rad/sec]

    Raises
    ------
    IOError
        Path can not be opened for reading
    AssertionError
        Parsing errors when reading

    """
    inp = open(file_swd, 'rb')

    def get_int():
        return unpack('<i', inp.read(4))[0]

    def get_float():
        return unpack('<f', inp.read(4))[0]

    def get_str(nstr):
        return unpack('<{0}s'.format(nstr), inp.read(nstr))[0]

    magic_number = get_float()
    assert abs(magic_number - 37.0221) < 0.0001, ("Wrong magic number in %" % file_swd)

    fmt = get_int()
    assert fmt == 100, ("Unknown format specifier: %i in file %s" % (fmt, file_swd))

    shp = get_int()
    assert shp == 6, ("Only shape 6 is supported in this function. (actual shp=%i)" % shp)

    amp = get_int()
    prog = get_str(30)
    date = get_str(20)
    nid = get_int()
    cid = get_str(nid)
    grav = get_float()
    Lscale = get_float()
    nstrip = get_int()
    nsteps = get_int()
    dt = get_float()
    order = get_int()
    nwaves = get_int()
    depth = get_float()

    amps = np.empty(nwaves)
    kwaves = np.empty(nwaves)
    omegas = np.empty(nwaves)
    Twaves = np.empty(nwaves)
    Lwaves = np.empty(nwaves)
    dirs_deg = np.empty(nwaves)
    phases_deg = np.empty(nwaves)

    for i in range(nwaves):

        amps[i] = get_float()
        kwaves[i] = get_float()
        dirs_deg[i] = get_float() * 180.0 / np.pi
        phases_deg[i] = get_float() * 180.0 / np.pi

        if depth < 0.0:
            omegas[i] = np.sqrt(kwaves[i] * grav)
        else:
            omegas[i] = np.sqrt(kwaves[i] * grav * np.tanh(kwaves[i] * depth))

        Twaves[i] = 2.0 * np.pi / omegas[i]
        Lwaves[i] = 2.0 * np.pi / kwaves[i]

    inp.close()

    res = {"n" : nwaves,
           "depth" : depth,
           "grav" : grav,
           "amps" : amps,
           "kwaves" : kwaves,
           "omegas" : omegas,
           "Twaves" : Twaves,
           "Lwaves" : Lwaves,
           "dirs_deg" : dirs_deg,
           "phases_deg" : phases_deg}

    return res


def omega2kwave(omega, depth, grav=9.81):
    """
    Solve the linear dispersion relation close to machine precision::

        omega**2 = kwave * grav * tanh(kwave*depth)

    Parameters
    ----------
    omega : float
        Wave oscillation frequency [rad/s]
    depth : float
        Constant water depth. [m] (<0 indicates infinite depth)
    grav : float, optional
        Acceleration of gravity [m/s^2]

    Returns
    -------
    float
        Wave number (kwave) [1/m]

    Raises
    ------
    None

    """

    if depth < 0.0:
        return omega**2 / grav

    # Solve equivalent equation system: c == y * tanh(y),  kwave = y / depth
    c = depth * omega**2 / grav

    # High accuracy fix point schemes
    if c > 2.5:
        def f(y0):
            # y == c/tanh(y)
            # tanh(y) = 1 - eps, Near y=y0 the RHS is almost c.
            # Solve y== c / tanh(y0)
            return c / np.tanh(y0)
    else:
        def f(y0):
            # y*(k0 + k1*(y-y0)) == c*(k0 + k1*(y-y0))/tanh(y0)
            # tanh(y0) = k0 + k1*(y-y0) + ...
            # Near y=y0 the RHS is almost c.
            # Solve y*(k0 + k1*(y-y0)) == c for y
            k0 = np.tanh(y0)
            k1 = 1.0 - k0**2
            b = k0 - k1 * y0
            return 0.5 * (-b + np.sqrt(b**2 + 4.0 * c * k1)) / k1

    # Fist initial guess (MIT lecture notes... 4 digits accuracy)
    if c > 2.4:
        # Deeper water...
        y = c * (1.0 + 3.0 * np.exp(-2 * c) - 12.0 * np.exp(-4 * c))
        # using fixed point iteration: y <= c + y - y * tanh(y)
    else:
        # Shallower water...
        y = np.sqrt(c) * (1.0 + 0.169 * c + 0.031 * c ** 2)
        # using fixed point iteration: y <= sqrt(c * y / tanh(y))

    y_prev = -1.0
    while abs(y - y_prev) > 100 * np.finfo(y).eps:
        y_prev = y
        y = f(y)
    kwave = y / depth

    return kwave
