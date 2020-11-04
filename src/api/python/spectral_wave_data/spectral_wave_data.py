# -*- coding: utf-8 -*-
    
"""
:platform: Linux, Windows, python 2.7 and 3.x
:synopsis: Python wrapper for the spectral_wave_data ocean wave model.

Author  - Jens Bloch Helmers, DNVGL
Created - 2019-08-11
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from ctypes import POINTER, c_double, cast

__all__ = ['SpectralWaveData', 'SwdError', 'SwdFileCantOpenError',
           'SwdFileBinaryError', 'SwdFileDataError', 'SwdInputValueError',
           'SwdAllocateError']

from .swd_c_interface import swdlib, vecswd, vecphi2ndswd, vecelev2ndswd


class SwdError(Exception):
    pass


class SwdFileCantOpenError(SwdError):
    pass


class SwdFileBinaryError(SwdError):
    pass


class SwdFileDataError(SwdError):
    pass


class SwdInputValueError(SwdError):
    pass


class SwdAllocateError(SwdError):
    pass


class SpectralWaveData(object):
    """An object of SpectraWaveData is an instance of the SWD ocean wave model.

    Raises
    ------
    SwdError
         Base class for the following exceptions
    SwdFileCantOpenError
         Typical if SWD file is not an existing file
    SwdFileBinaryError
         SWD file does not apply float/little-endian
    SwdFileDataError
         Error during reading and checking data
    SwdInputValueError
         Input arguments for class methods are not sound
    SwdAllocateError
         Not able to allocate internal SWD storage

    Attributes
    ----------
    No public attributes

    """

    def __init__(self, file_swd, x0=0.0, y0=0.0, t0=0.0, beta=0.0, rho=1025.0, nsumx=-1, 
                 nsumy=-1, impl=0, ipol=0, norder=0, dc_bias=False):
        """Constructor

        Parameters
        ----------
        file_swd : str
            The name of the swd file defining the ocean waves.
        x0, y0 : float, optional
            The origin of the application wave coordinate system relative to the
            SWD coordinate system. [m]
        t0 : float, optional
            The SWD time corresponding to t=0 in the application simulation. [s]
        beta : float, optional
            Rotation of the SWD x-axis relative to the application x-axis. [deg]
        rho : float, optional
            Density of water. [kg/m^3] (Only relevant for pressure calculations)
        nsumx, nsumy : int, optional
            Number of applied spectral components in x and y-directions. Higher frequency
            components in the swd file will not be applied. nsumx < 0 and nsumy < 0
            indicates that all components from the swd file will be applied.
        impl : int, optional
            Index to request actual derived SWD class
                *  0 = Apply a default class based on the content of this SWD-file.
                * <0 = Apply an in-house and experimental implementation
                * >0 = Apply a specific implementations from the open software
        ipol : int, optional
            Index to request actual temporal interpolation scheme
                *  0 = Default (C^2 continous scheme)
                *  1 = C^1 continous scheme
                *  2 = C^3 continous scheme
        norder : int, optional
             Expansion order to apply in kinematics for z>0
                *  0 = Apply expansion order specified on swd file (default)
                * <0 = Apply exp(kj z)
                * >0 = Apply expansion order = norder
        dc_bias : bool, optional
             Control application of zero-frequency bias present in SWD file
                * False = Suppress contribution from zero frequency amplitudes (default)
                * True  = Apply zero frequency amplitudes from SWD file.

        Raises
        ------
        SwdError
             Base class for the following exceptions
        SwdFileCantOpenError
             Typical if SWD file is not an existing file
        SwdFileBinaryError
             SWD file does not apply float/little-endian
        SwdFileDataError
             Error during reading and checking data
        SwdInputValueError
             Input arguments for class methods are not sound
        SwdAllocateError
             Not able to allocate internal SWD storage

        Examples
        --------
        >>> from spectral_wave_data import SpectralWaveData
        >>> swd = SpectralWaveData('my_waves.swd', x0=0.0, y0=0.0, t0=0.0, beta=180.0)

        """
        self.obj = swdlib.swd_api_allocate(file_swd.encode('ascii'), x0, y0, t0, beta,
                                           rho, nsumx, nsumy, impl, ipol, norder, dc_bias)
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            if id == 1001:
                raise SwdFileCantOpenError(msg)
            elif id == 1002:
                raise SwdFileBinaryError(msg)
            elif id == 1003:
                raise SwdFileDataError(msg)
            elif id == 1004:
                raise SwdInputValueError(msg)
            elif id == 1005:
                raise SwdAllocateError(msg)
            else:
                raise SwdError(msg)
        self._alive = True

    def update_time(self, time):
        """Set the current time for all kinematic calculations.

        .. note::

           This method must be called at least once before any kinematic
           calculations can be done.

        Parameters
        ----------
        time : float
            Time as defined in the application program [s]

        Returns
        -------
        None

        Raises
        ------
        SwdError
             Base class for the following exceptions
        SwdFileDataError
             Error during reading the SWD file
        SwdInputValueError
             Input arguments are not sound

        Examples
        --------
        >>> swd.update_time(time=73.241)    # Time as defined in the application

        """
        swdlib.swd_api_update_time(self.obj, time)
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1003:
                raise SwdFileDataError(msg)
            elif id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)

    def phi(self, x, y, z):
        """Calculates the velocity potential at the actual location. It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        float
            Velocity potential at (x,y,z) [m^2/s]

        Raises
        ------
        None

        Examples
        --------
        >>> print("potential at (x,y,z) = ", swd.phi(x,y,z))

        """
        res = swdlib.swd_api_phi(self.obj, x, y, z)
        return res

    def stream(self, x, y, z):
        """Calculates the stream function at the actual location. It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        float
            Stream function at (x,y,z)

        Raises
        ------
        None

        Examples
        --------
        >>> print("Stream function at (x,y,z) = ", swd.stream(x,y,z))

        """
        res = swdlib.swd_api_stream(self.obj, x, y, z)
        return res

    def phi_t(self, x, y, z):
        """Calculates the time derivative of the wave potential at the actual
        location (Euler derivative). It is assumed that the current time has
        been set using the method :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        float
            Time derivative (Euler) of wave potential at (x,y,z) [m^2/s^2]

        Raises
        ------
        None

        Examples
        --------
        >>> print("delta(phi)/delta(t) at (x,y,z) = ", swd.phi_t(x,y,z))

        """
        res = swdlib.swd_api_phi_t(self.obj, x, y, z)
        return res

    def grad_phi(self, x, y, z):
        """Calculates the particle velocity. The velocity components are
        with respect to the application coordinate system. It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        vecswd (struct with x, y and z float attributes)
            Particle velocity at (x,y,z) [m/s]

        Raises
        ------
        None

        Examples
        --------
        >>> vel = swd.grad_phi(x, y, z)
        >>> print("velocity in x-dir = ", vel.x)
        >>> print("velocity in y-dir = ", vel.y)
        >>> print("velocity in z-dir = ", vel.z)

        """
        res = swdlib.swd_api_grad_phi(self.obj, x, y, z)
        return res      # res.x, res.y, res.z

    def grad_phi_2nd(self, x, y, z):
        """Calculates the 2nd order derivatives of the wave potential.
        Gradients are with with respect to the application coordinate system.
        It is assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        vecphi2ndswd (struct with xx, xy, xz, yy, yz, zz float attributes)
            2nd order derivatives of velocity potential [1/s]

        Raises
        ------
        None

        Examples
        --------
        >>> res = swd.grad_phi_2nd(x, y, z)
        >>> print("phi_xx = ", res.xx)
        >>> print("phi_xy = ", res.xy)
        >>> print("phi_xz = ", res.xz)
        >>> print("phi_yy = ", res.yy)
        >>> print("phi_yz = ", res.yz)
        >>> print("phi_zz = ", res.zz)

        """
        res = swdlib.swd_api_grad_phi_2nd(self.obj, x, y, z)
        return res       # res.xx, res.xy, res.xz, res.yy, res.yz, res.zz

    def acc_euler(self, x, y, z):
        """Calculates the Euler acceleration (grad(phi_t)) at (x,y,z).
        Components are with respect to the application coordinate system.
        It is assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        vecswd (struct with x, y and z float attributes)
            Euler acceleration at (x,y,z) [m/s^2]

        Raises
        ------
        None

        Examples
        --------
        >>> acc = swd.acc_euler(x, y, z)
        >>> print("Euler acceleration in x-dir = ", acc.x)
        >>> print("Euler acceleration in y-dir = ", acc.y)
        >>> print("Euler acceleration in z-dir = ", acc.z)

        """
        res = swdlib.swd_api_acc_euler(self.obj, x, y, z)
        return res      # res.x, res.y, res.z

    def acc_particle(self, x, y, z):
        """Calculates the particle acceleration at (x,y,z).
        Components are with respect to the application coordinate system.
        It is assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        vecswd (struct with x, y and z float attributes)
            Particle acceleration at (x,y,z) [m/s^2]

        Raises
        ------
        None

        Examples
        --------
        >>> acc = swd.acc_particle(x, y, z)
        >>> print("Particle acceleration in x-dir = ", acc.x)
        >>> print("Particle acceleration in y-dir = ", acc.y)
        >>> print("Particle acceleration in z-dir = ", acc.z)

        """
        res = swdlib.swd_api_acc_particle(self.obj, x, y, z)
        return res      # res.x, res.y, res.z

    def elev(self, x, y):
        """Calculates the wave elevation at (x,y). It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y : float
            Horizontal position as defined in the application program [m]

        Returns
        -------
        float
            Wave elevation at (x,y) [m]

        Raises
        ------
        None

        Examples
        --------
        >>> print("Wave elevation at (x,y) = ", swd.elev(x,y))

        """
        res = swdlib.swd_api_elev(self.obj, x, y)
        return res

    def elev_t(self, x, y):
        """Calculates the time derivative of the wave elevation
        at earth fixed location. It is assumed that the current
        time has been set using the method :meth:`update_time`.

        Parameters
        ----------
        x, y : float
            Horizontal position as defined in the application program [m]

        Returns
        -------
        float
            Time derivative of wave elevation at (x,y)  [m/s]

        Raises
        ------
        None

        Examples
        --------
        >>> print("Time derivative of elevation at (x,y) = ", swd.elev_t(x,y))

        """
        res = swdlib.swd_api_elev_t(self.obj, x, y)
        return res

    def grad_elev(self, x, y):
        """Calculates the gradients of the wave elevation. The gradients
        are with respect to the application coordinate system. It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y : float
            Horizontal location as defined in the application program [m]

        Returns
        -------
        vecswd (struct with x, y and z float attributes)
            Spatial gradients of surface elevation at (x,y) [-]

        Raises
        ------
        None

        Examples
        --------
        >>> res = swd.grad_elev(x, y)
        >>> print("d/dx(elevation) = ", res.x)
        >>> print("d/dy(elevation) = ", res.y)
        >>> res.z
        0.0

        """
        res = swdlib.swd_api_grad_elev(self.obj, x, y)
        return res      # res.x, res.y, res.z=0

    def grad_elev_2nd(self, x, y):
        """Calculates the 2nd order gradients of the wave elevations. The gradients
        are with respect to the application coordinate system. It is
        assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y : float
            Horizontal location as defined in the application program [m]

        Returns
        -------
        vecelev2ndswd (struct with xx, xy and yy float attributes)
            2nd order derivatives of wave elevation at (x,y) [1/m]

        Raises
        ------
        None

        Examples
        --------
        >>> res = swd.grad_elev_2nd(x, y)
        >>> print("elvation_xx = ", res.xx)
        >>> print("elvation_xy = ", res.xy)
        >>> print("elvation_yy = ", res.yy)

        """
        res = swdlib.swd_api_grad_elev_2nd(self.obj, x, y)
        return res       # res.xx, res.xy, res.yy

    def bathymetry(self, x, y):
        """Calculates the vertical distance from z=0 to the sea floor at (x, y).

        Parameters
        ----------
        x, y : float
            Horizontal position as defined in the application program [m]

        Returns
        -------
        float
            Local water depth at (x,y)  (<0 indicates infinite depth) [m]

        Raises
        ------
        None

        Examples
        --------
        >>> print("Local water depth at (x,y) = ", swd.bathymetry(x,y))

        """
        res = swdlib.swd_api_bathymetry(self.obj, x, y)
        return res

    def bathymetry_nvec(self, x, y):
        """Return the unit normal vector of the sea floor at (x,y).
        The orientation of the vector is into the ocean and with
        respect to the application coordinate system.

        Parameters
        ----------
        x, y : float
            Horizontal location as defined in the application program [m]

        Returns
        -------
        vecswd (struct with x, y and z float attributes)
            Unit normal vector on sea floor at (x,y)

        Raises
        ------
        None

        Examples
        --------
        >>> nvec = swd.bathymetry_nvec(x, y)
        >>> print("n_x = ", nvec.x)
        >>> print("n_y = ", nvec.y)
        >>> print("n_z = ", nvec.z)  # Typical close to 1

        """
        res = swdlib.swd_api_bathymetry_nvec(self.obj, x, y)
        return res    # res.x, res.y, res.z

    def pressure(self, x, y, z):
        """Calculates the complete nonlinear Bernoulli-pressure (Const=0).
        It is assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]

        Returns
        -------
        float
            Pressure at (x,y,z)  [Pa]

        Raises
        ------
        None

        Examples
        --------
        >>> print("Total pressure at (x,y,z) = ", swd.pressure(x,y,z))

        """
        res = swdlib.swd_api_pressure(self.obj, x, y, z)
        return res

    def convergence(self, x, y, z, csv):
        """For a specific location create a CSV-file on how velocity, elevation
        and pressure converge as a function of number of spectral components.
        It is assumed that the current time has been set using the method
        :meth:`update_time`.

        Parameters
        ----------
        x, y, z : float
            Position as defined in the application program [m]
        csv : str
            Name of requested CSV-file to be generated.

        Returns
        -------
        None

        Raises
        ------
        SwdError
            Base class for the following exceptions
        SwdFileCantOpenError
            If an existing CSV-file with the same name is locked.

        Examples
        --------
        >>> swd.convergence(x, y, z, 'convergence_data_at_xyz.csv')

        """
        swdlib.swd_api_convergence(self.obj, x, y, z, csv.encode('ascii'))
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1001:
                raise SwdFileCantOpenError(msg)
            else:
                raise SwdError(msg)

    def strip(self, tmin, tmax, file_swd):
        """Create a new SWD file containing only the time steps within the
        time window [tmin, tmax].

        Parameters
        ----------
        tmin, tmax : float
            Time window as defined in the application program [s]
        file_swd : str
            Name of new SWD file containing only the actual time window

        Returns
        -------
        None

        Raises
        ------
        SwdError
            Base class for the following exceptions
        SwdFileCantOpenError
            If not able to open new SWD file
        SwdFileDataError
            Error during reading data from existing SWD file
        SwdInputValueError
            The time window is outside the content of existing SWD file.

        Examples
        --------
        >>> swd.strip(tmin=850.0, tmax=950.0, file_swd='freak_wave_at_850_950.swd')

        """
        swdlib.swd_api_strip(self.obj, tmin, tmax, file_swd.encode('ascii'))
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1001:
                raise SwdFileCantOpenError(msg)
            elif id == 1003:
                raise SwdFileDataError(msg)
            elif id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)

    def __getitem__(self, key):
        """self[key] is an alias for calling the method :meth:`get`

        Examples
        --------
        >>> swd = SpectralWaveData('my_waves.swd', ...)
        >>> swd['tmax']     # Return max allowed application time

        """
        return self.get(key)

    def get(self, key):
        """Extract metadata from the SWD object

        Parameters
        ----------
        key : str
            A key to identify requested metadata. A key is either:
                 * A relevant parameter from the SWD-file format description.
                 * A constructor parameter.
                 * A key from the table below.

            ::

               Key         Returned Value
               ----------------------------------------------------------------------
               'version'   The repository version number of spectral-wave-data
                           applied in this Python distribution.
               'class'     Name of the specialized SWD-class handling this object.
               'tmax'      Maximum allowed application time consistent with the
                           header of the SWD file.
               'lmin'      Shortest wave length component. [m]
               'lmax'      Longest wave length component. [m]
               'sizex'     Periodic length of wave domain in x-direction (swd-system).
               'sizey'     Periodic length of wave domain in y-direction (swd-system).

        Returns
        -------
        float, int, bool, or str
            Value of actual metadata

        Raises
        ------
        SwdInputValueError
            Key does not correspond to relevant metadata for this object.

        Examples
        --------
        >>> print("swd version number = ", swd.get('version'))
        >>> print("Max allowed user time = ", swd.get('tmax'))
        >>> print("Size of periodic domain in x-dir (SWD) = ", swd.get('sizex'))
        >>> print("Actual spectral shape model = ", swd.get('shp'))
        >>> print("Applied SWD-class = ", swd.get('class'))

        """
        key_c = key.encode('ascii')

        if key in ['file', 'file_swd', 'version', 'class', 'cid', 'prog', 'date']:
            res = swdlib.swd_api_get_chr(self.obj, key_c).decode()
        elif key in ['magic', 'grav', 'lscale', 'dt', 'dk', 'dkx', 'dky',
                      'd', 'depth', 'tmax', 'lmin', 'lmax', 
                      'sizex', 'sizey', 't0', 'x0', 'y0', 'beta', 'rho']:
            res = swdlib.swd_api_get_real(self.obj, key_c)
        elif key in ['dc_bias']:
            res = swdlib.swd_api_get_bool(self.obj, key_c)
        else:
            res = swdlib.swd_api_get_int(self.obj, key_c)

        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)  # Just in case....

        return res

    def close(self):
        """Manual destructor closing the object and the SWD-file.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        None

        """
        if self._alive is True:
            swdlib.swd_api_close(self.obj)
            self._alive = False

    def elev_fft(self, nx_fft=-1, ny_fft=-1):
        """Calculates the wave elevation at a regular FFT-grid with
        dimensions (nx_fft, ny_fft). It is assumed that the current time has 
        been set using the method :meth:`update_time`.

        Parameters
        ----------
        nx_fft, ny_fft : int, optional
            Output dimensions.

        Returns
        -------
        ndarray
            Wave elevation on a regular grid.

        Raises
        ------
        SwdInputValueError
            Invalid grid dimensions (nx_fft, ny_fft).

        Examples
        --------
        

        """
        # get the fortran-array-object (see ISO_Fortran_binding.h/ISO_Fortran_binding.py)
        CFI_obj = swdlib.swd_api_elev_fft(self.obj, nx_fft, ny_fft)    
        
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)

        # array dimensions
        nx_out = CFI_obj.contents.dim[0].extent
        ny_out = CFI_obj.contents.dim[1].extent

        data_pointer = cast(CFI_obj.contents.base_addr, POINTER(c_double))
        res = np.ctypeslib.as_array(data_pointer, shape=(ny_out, nx_out)).copy()

        swdlib.swd_api_fft_deallocate(CFI_obj)
        
        return res.T

    def x_fft(self, nx_fft=-1):
        """Returns the x-grid corresponding to the FFT-based routines *_fft
        called with the same nx_fft.

        Parameters
        ----------
        nx_fft : int, optional
            Grid dimensions.

        Returns
        -------
        ndarray
            x-values on the FFT-grid.

        Raises
        ------
        SwdInputValueError
            Invalid grid dimension nx_fft.

        Examples
        --------
        

        """
        # get the fortran-array-object (see ISO_Fortran_binding.h/ISO_Fortran_binding.py)
        CFI_obj = swdlib.swd_api_x_fft(self.obj, nx_fft)    
        
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)

        # array dimensions
        nx_out = CFI_obj.contents.dim[0].extent

        data_pointer = cast(CFI_obj.contents.base_addr, POINTER(c_double))
        res = np.ctypeslib.as_array(data_pointer, shape=(nx_out, )).copy()

        swdlib.swd_api_fft_deallocate(CFI_obj)
        
        return res
    
    def y_fft(self, ny_fft=-1):
        """Returns the y-grid corresponding to the FFT-based routines *_fft
        called with the same ny_fft.

        Parameters
        ----------
        ny_fft : int, optional
            Grid dimensions.

        Returns
        -------
        ndarray
            y-values on the FFT-grid.

        Raises
        ------
        SwdInputValueError
            Invalid grid dimension ny_fft.

        Examples
        --------
        

        """
        # get the fortran-array-object (see ISO_Fortran_binding.h/ISO_Fortran_binding.py)
        CFI_obj = swdlib.swd_api_y_fft(self.obj, ny_fft)    
        
        if swdlib.swd_api_error_raised(self.obj):
            id = swdlib.swd_api_error_get_id(self.obj)
            msg = swdlib.swd_api_error_get_msg(self.obj).decode()
            swdlib.swd_api_error_clear(self.obj) # To simplify safe recovery...
            if id == 1004:
                raise SwdInputValueError(msg)
            else:
                raise SwdError(msg)

        # array dimensions
        ny_out = CFI_obj.contents.dim[0].extent

        data_pointer = cast(CFI_obj.contents.base_addr, POINTER(c_double))
        res = np.ctypeslib.as_array(data_pointer, shape=(ny_out, )).copy()

        swdlib.swd_api_fft_deallocate(CFI_obj)
        
        return res
