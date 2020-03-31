from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import os
from struct import pack

import sympy as sp
import numpy as np

from spectral_wave_data.tools import airy

from symbols import x, y, z, t

class Shape6:

    def __init__(self, amps, dirs, phases, kwaves, depth, order, grav, sys):
        assert order < 3
        assert len(amps) == len(dirs) == len(phases) == len(kwaves)
        self.n = len(amps)
        self.amps = np.array(amps)
        self.dirs = np.array(dirs)  # Radians
        self.phases = np.array(phases) # Radians
        self.kwaves = np.array(kwaves)
        self.kwx = self.kwaves * np.cos(self.dirs)
        self.kwy = self.kwaves * np.sin(self.dirs)
        self.long_crested = np.allclose(self.dirs, self.dirs[0])
        self.grav = grav
        self.depth = depth
        self.order = order
        if depth < 0.0:
            self.omegas = np.sqrt(self.kwaves * self.grav)
        else:
            self.omegas = np.sqrt(self.kwaves * self.grav * np.tanh(self.kwaves * depth))
        self.sys = sys
        # Express swd coordinates using application coordinates (x,y,z,t)
        self.xswd, self.yswd, self.zswd, self.tswd = sys.app2swd()

        fun_phi = 0
        fun_varphi = 0
        fun_elev = 0
        for j in range(self.n):
            arg = self.omegas[j] * self.tswd - \
                  self.kwx[j] * self.xswd - \
                  self.kwy[j] * self.yswd + self.phases[j]
            phi_amp = -self.grav * self.amps[j] / self.omegas[j]
            fun_phi += phi_amp * self.Zfun(j) * sp.sin(arg)
            fun_varphi += -phi_amp * self.Zhfun(j) * sp.cos(arg)
            fun_elev += self.amps[j] * sp.cos(arg)

        self.f_phi = fun_phi
        if self.long_crested:
            self.f_stream = fun_varphi
        else:
            self.f_stream = sp.Integer(0)

        self.f_phi_t = sp.diff(fun_phi, t)
        self.f_phi_x = sp.diff(fun_phi, x)
        self.f_phi_y = sp.diff(fun_phi, y)
        self.f_phi_z = sp.diff(fun_phi, z)

        self.f_phi_tx = sp.diff(fun_phi, t, x)
        self.f_phi_ty = sp.diff(fun_phi, t, y)
        self.f_phi_tz = sp.diff(fun_phi, t, z)

        self.f_phi_xx = sp.diff(fun_phi, x, x)
        self.f_phi_xy = sp.diff(fun_phi, x, y)
        self.f_phi_xz = sp.diff(fun_phi, x, z)
        self.f_phi_yy = sp.diff(fun_phi, y, y)
        self.f_phi_yz = sp.diff(fun_phi, y, z)
        self.f_phi_zz = sp.diff(fun_phi, z, z)

        # Linearize above functions with respect to z=0
        self.f_phi_lin = self.f_phi.series(z, 0, 2).removeO()
        self.f_stream_lin = self.f_stream.series(z, 0, 2).removeO()
        self.f_phi_t_lin = self.f_phi_t.series(z, 0, 2).removeO()
        self.f_phi_x_lin = self.f_phi_x.series(z, 0, 2).removeO()
        self.f_phi_y_lin = self.f_phi_y.series(z, 0, 2).removeO()
        self.f_phi_z_lin = self.f_phi_z.series(z, 0, 2).removeO()
        self.f_phi_tx_lin = self.f_phi_tx.series(z, 0, 2).removeO()
        self.f_phi_ty_lin = self.f_phi_ty.series(z, 0, 2).removeO()
        self.f_phi_tz_lin = self.f_phi_tz.series(z, 0, 2).removeO()
        self.f_phi_xx_lin = self.f_phi_xx.series(z, 0, 2).removeO()
        self.f_phi_xy_lin = self.f_phi_xy.series(z, 0, 2).removeO()
        self.f_phi_xz_lin = self.f_phi_xz.series(z, 0, 2).removeO()
        self.f_phi_yy_lin = self.f_phi_yy.series(z, 0, 2).removeO()
        self.f_phi_yz_lin = self.f_phi_yz.series(z, 0, 2).removeO()
        self.f_phi_zz_lin = self.f_phi_zz.series(z, 0, 2).removeO()

        self.f_elev = fun_elev
        self.f_elev_t = sp.diff(fun_elev, t)

        self.f_elev_x = sp.diff(fun_elev, x)
        self.f_elev_y = sp.diff(fun_elev, y)

        self.f_elev_xx = sp.diff(fun_elev, x, x)
        self.f_elev_xy = sp.diff(fun_elev, x, y)
        self.f_elev_yy = sp.diff(fun_elev, y, y)

    def z_apply(self, x_app, y_app, z_app, t_app):
        if self.order < 0:
            z = z_app
        elif self.order == 0:
            z = min(z_app, 0)
        elif self.order == 1:
            z = z_app
        elif self.order == 2:
            elev = self.elev(x_app, y_app, t_app)
            z = z_app - elev
            if self.depth > 0:
                z = z/(1 + elev/self.depth)
        return z

    def Zfun(self, j):
        if self.depth < 0:
            return sp.exp(self.kwaves[j] * self.zswd)
        else:
            return sp.cosh(self.kwaves[j] * (self.zswd + self.depth)) / \
                   sp.cosh(self.kwaves[j] * self.depth)

    def Zhfun(self, j):
        if self.depth < 0:
            return sp.exp(self.kwaves[j] * self.zswd)
        else:
            return sp.sinh(self.kwaves[j] * (self.zswd + self.depth)) / \
                   sp.cosh(self.kwaves[j] * self.depth)

    def bathymetry(self, x_app, y_app):
        return self.depth

    def bathymetry_nvec(self, x_app, y_app):
        return {'x': 0.0, 'y': 0.0, 'z': 1.0}

    def phi(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            return sp.N(self.f_phi_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
        else:
            return sp.N(self.f_phi.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))

    def stream(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            return sp.N(self.f_stream_lin.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
        else:
            return sp.N(self.f_stream.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))

    def phi_t(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            return sp.N(self.f_phi_t_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
        else:
            return sp.N(self.f_phi_t.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))

    def grad_phi(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            phi_x = sp.N(self.f_phi_x_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_y = sp.N(self.f_phi_y_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_z = sp.N(self.f_phi_z_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
        else:
            phi_x = sp.N(self.f_phi_x.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_y = sp.N(self.f_phi_y.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_z = sp.N(self.f_phi_z.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
        return {'x':phi_x, 'y':phi_y, 'z':phi_z}

    def grad_phi_2nd(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            phi_xx = sp.N(self.f_phi_xx_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_xy = sp.N(self.f_phi_xy_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_xz = sp.N(self.f_phi_xz_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_yy = sp.N(self.f_phi_yy_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_yz = sp.N(self.f_phi_yz_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_zz = sp.N(self.f_phi_zz_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
        else:
            phi_xx = sp.N(self.f_phi_xx.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_xy = sp.N(self.f_phi_xy.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_xz = sp.N(self.f_phi_xz.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_yy = sp.N(self.f_phi_yy.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_yz = sp.N(self.f_phi_yz.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_zz = sp.N(self.f_phi_zz.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
        return {'xx':phi_xx, 'xy':phi_xy, 'xz':phi_xz,
                'yy':phi_yy, 'yz':phi_yz, 'zz':phi_zz}

    def acc_euler(self, x_app, y_app, z_app, t_app):
        z_eval = self.z_apply(x_app, y_app, z_app, t_app)
        if self.order == 1 and z_eval > 0:
            phi_tx = sp.N(self.f_phi_tx_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_ty = sp.N(self.f_phi_ty_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
            phi_tz = sp.N(self.f_phi_tz_lin.subs({x:x_app, y:y_app, z:z_eval, t:t_app}))
        else:
            phi_tx = sp.N(self.f_phi_tx.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_ty = sp.N(self.f_phi_ty.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
            phi_tz = sp.N(self.f_phi_tz.subs({x: x_app, y: y_app, z: z_eval, t: t_app}))
        return {'x':phi_tx, 'y':phi_ty, 'z':phi_tz}

    def acc_particle(self, x_app, y_app, z_app, t_app):
        vel = self.grad_phi(x_app, y_app, z_app, t_app)
        aeul = self.acc_euler(x_app, y_app, z_app, t_app)
        g2nd = self.grad_phi_2nd(x_app, y_app, z_app, t_app)
        ax = sp.N(aeul['x'] + vel['x'] * g2nd['xx'] + vel['y'] * g2nd['xy'] + vel['z'] * g2nd['xz'])
        ay = sp.N(aeul['y'] + vel['x'] * g2nd['xy'] + vel['y'] * g2nd['yy'] + vel['z'] * g2nd['yz'])
        az = sp.N(aeul['z'] + vel['x'] * g2nd['xz'] + vel['y'] * g2nd['yz'] + vel['z'] * g2nd['zz'])
        return {'x':ax, 'y':ay, 'z':az}

    def pressure(self, x_app, y_app, z_app, t_app, rho, grav):
        vel = self.grad_phi(x_app, y_app, z_app, t_app)
        phit = self.phi_t(x_app, y_app, z_app, t_app)
        # Note that the hydrostatic component is always evaluated at physical location:
        return sp.N(-rho * phit - rho * grav * z_app - 0.5 * rho * (vel['x']**2 + vel['y']**2 + vel['z']**2))

    def elev(self, x_app, y_app, t_app):
        return sp.N(self.f_elev.subs({x:x_app, y:y_app, t:t_app}))

    def elev_t(self, x_app, y_app, t_app):
        return sp.N(self.f_elev_t.subs({x:x_app, y:y_app, t:t_app}))

    def grad_elev(self, x_app, y_app, t_app):
        elv_x = sp.N(self.f_elev_x.subs({x:x_app, y:y_app, t:t_app}))
        elv_y = sp.N(self.f_elev_y.subs({x:x_app, y:y_app, t:t_app}))
        elv_z = 0
        return {'x':elv_x, 'y':elv_y, 'z':elv_z}

    def grad_elev_2nd(self, x_app, y_app, t_app):
        elv_xx = sp.N(self.f_elev_xx.subs({x:x_app, y:y_app, t:t_app}))
        elv_xy = sp.N(self.f_elev_xy.subs({x:x_app, y:y_app, t:t_app}))
        elv_yy = sp.N(self.f_elev_yy.subs({x:x_app, y:y_app, t:t_app}))
        return {'xx':elv_xx, 'xy':elv_xy, 'yy':elv_yy}


    def write_swd(self, file_swd):
        airy.write_swd(file_swd, self.amps, self.dirs, self.phases, kwaves=self.kwaves,
                       depth=self.depth, grav=self.grav, is_deg_dirs=False, is_deg_phases=False)
