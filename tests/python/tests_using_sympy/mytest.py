from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import symbols

import tfun
from corsys import CorSys
from shape_1 import Shape1


# Temporal interpolation scheme as applied in official software
ipol = 0

# Number of non-bias spectral components
n_swd = 3

# Seed for random number generator
seed = 0

# Definition of application coordinate system
sys = CorSys(x0=0, y0=0, t0=0, beta=0)

# Define all temporal functions using polynomial of proper order
# related to the the ipol-scheme applied in the numerical method.
# Consequently, all temporal calculations should be exact.
# The polynomial coefficients are random in [-1.0, 1.0]
cfuns, hfuns = tfun.create_Tfuns_1d(n_swd, ipol=ipol, seed=seed)

# Define the fields...
shp = Shape1(dk=0.1, n=n_swd, cfuns=cfuns, hfuns=hfuns, order=2, sys=sys)

# The explicit analytical field functions with respect to the application
# system can be extracted in this way...
print("f_phi = ", shp.f_phi)
print("f_phi_t = ", shp.f_phi_t)

# Evaluations of symbolic expressions to be compared with official software
print('phi = ', shp.phi(x_app=0.3, y_app=-2.3, z_app=-0.2, t_app=0.2))
print('stream = ', shp.stream(x_app=0.3, y_app=-2.3, z_app=-0.2, t_app=0.2))
print('phi_t = ', shp.phi_t(x_app=0.3, y_app=-2.3, z_app=-0.2, t_app=0.2))
print('grad(phi) = ', shp.grad_phi(x_app=0.3, y_app=-2.3, z_app=-0.2, t_app=0.2))
print('elev = ', shp.elev(x_app=0.3, y_app=-2.3, t_app=0.2))

# Dump fields to swd file
shp.write_swd('test.swd', dt=0.1, nsteps=21)