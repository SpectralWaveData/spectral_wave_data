from __future__ import print_function   # Supports print() in Python 2.7 too

from spectral_wave_data import SpectralWaveData, SwdError  # Supports Python 2.x and 3.x


# First we define some constructor parameters...

# Relate the application and SWD coordinate systems
x0 = 0.0; y0 = 0.0; t0 = 0.0; beta = 0.0

# Optional: impl=0 means apply recommended implementation based on header of swd file.
impl = 0

# Optional: Density of water needed for pressure calculations
rho = 1025.0

# Optional: number of spectral components to apply in calculations.
nsumx = -1; nsumy = -1   # Negative: apply all components from swd-file.

# Optional: Select temporal interpolation scheme
ipol = 0   # C^2 continous

# Optional: Select expansion order for calculating kinematics above z=0
norder = 0   # Apply same order as specified in SWD file

# Optional: Control handling of zero-frequency components in SWD file
dc_bias = False   # Suppress contributions from zero-frequency

# Constructor
try:
    swd = SpectralWaveData('stokesTest_shallowWater.swd', x0, y0, t0, beta,
                           rho=rho, nsumx=nsumx, nsumy=nsumy, impl=impl,
                           ipol=ipol, norder=norder, dc_bias=dc_bias)
except SwdError as e:
    print(str(e))
    # Do whatever necessary...
    raise

# Time domain simulation...
t = 0.0;  dt = 0.1;  tmax = 1.0;
while t < tmax:

    # Tell the swd object current application time...
    swd.update_time(t)

    # Application coordinate where we need kinematics...
    x = 4.3;  y = 5.4;  z = -1.7;

    # Current set of provided kinematic quantities...

    # Velocity potential
    print("phi = ", swd.phi(x, y, z))

    # Stream function
    print("varphi = ", swd.stream(x, y, z))

    # partial(phi) / partial t
    print("phi_t = ", swd.phi_t(x, y, z))

    # particle velocity
    print("grad_phi = ", swd.grad_phi(x, y, z))

    # 2nd order gradients of phi
    print("phi_xx, phi_xy, phi_xz, phi_yy, phi_yz, phi_zz = ", swd.grad_phi_2nd(x, y, z))

    # Local (Euler) acceleration
    print("acc_euler = ", swd.acc_euler(x, y, z))

    # particle acceleration
    print("acc_particle = ", swd.acc_particle(x, y, z))

    # surface elevation
    print("zeta = ", swd.elev(x, y))

    # local time derivative of wave elevation
    print("\partial(zeta)/\partial t = ", swd.elev_t(x, y))

    # spatial gradient of wave elevation
    print("grad_zeta = ", swd.grad_elev(x, y))

    # 2nd order gradients of wave elevation
    print("zeta_xx, zeta_xy, zeta_yy = ", swd.grad_elev_2nd(x, y))

    # Total pressure
    print("pressure = ", swd.pressure(x, y, z))

    # Vertical distance from z=0 to sea floor (<0 if infinite)
    print("local depth = ", swd.bathymetry(x, y))

    # Unit normal vector of sea floor into the ocean at (x,y)
    print("nx, ny, nz = ", swd.bathymetry_nvec(x, y))

    # For a specific (x,y,z) at this t, return a CSV file on how particle velocity,
    # elevation and pressure converge as a function of number of spectral components
    swd.convergence(x, y, z, 'dump.csv')

    t += dt

# To save storage for an interesting event you may create a new SWD file
# containing only the time steps within a specified time window.
swd.strip(tmin=100.0, tmax=200.0, file_swd='my_new.swd')

# ===========================================================
# The meth swd.get(name) returns the value of parameter 'name'
# from the swd file. Three examples are given below...
# ===========================================================

# Extract the cid string from the SWD file
# (contains typical the content of the input file applied in the wave generator)
print("cid = ", swd.get('cid'))

# The shp parameter from the swd file
print("shp = ", swd.get('shp'))

# Time step in SWD file
print("dt = ", swd.get('dt'))

# Name of actual implementation class
print("cls_name = ", swd.get('class'))

# Close SWD file and free related memory
swd.close()
