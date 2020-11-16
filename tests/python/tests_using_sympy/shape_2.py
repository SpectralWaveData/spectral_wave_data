from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import subprocess
from struct import pack

import sympy as sp

from symbols import x, y, z, t


class Shape2:

    def __init__(self, dk, n, d, cfuns, hfuns, order, sys):
        assert len(cfuns) == len(hfuns) == n + 1
        assert order != 0
        self.dk = dk
        self.n = n
        self.d = d
        self.order = order
        self.sys = sys
        # Express swd coordinates using application coordinates (x,y,z,t)
        self.xswd, self.yswd, self.zswd, self.tswd = sys.app2swd()
        self.cfuns = cfuns
        self.hfuns = hfuns
        if self.order < 0:
            self.Zfun = self.Zfun_exact
            self.Zhfun = self.Zhfun_exact
        else:
            self.Zfun = self.Zfun_exact_taylor
            self.Zhfun = self.Zhfun_exact_taylor

        fun_phi = 0
        fun_varphi = 0
        fun_elev = 0
        for j in range(n + 1):
            cj = self.cfuns[j].fun(self.tswd)
            fun_phi += cj * self.Xfun(j) * self.Zfun(j)
            fun_varphi += cj * self.Xfun(j) * self.Zhfun(j)
            hj = self.hfuns[j].fun(self.tswd)
            fun_elev += hj * self.Xfun(j)

        self.f_phi = sp.re(fun_phi)
        self.f_stream = sp.im(fun_varphi)
        self.f_phi_t = sp.re(sp.diff(fun_phi, t))

        self.f_phi_x = sp.re(sp.diff(fun_phi, x))
        self.f_phi_y = sp.re(sp.diff(fun_phi, y))
        self.f_phi_z = sp.re(sp.diff(fun_phi, z))

        self.f_phi_tx = sp.re(sp.diff(fun_phi, t, x))
        self.f_phi_ty = sp.re(sp.diff(fun_phi, t, y))
        self.f_phi_tz = sp.re(sp.diff(fun_phi, t, z))

        self.f_phi_xx = sp.re(sp.diff(fun_phi, x, x))
        self.f_phi_xy = sp.re(sp.diff(fun_phi, x, y))
        self.f_phi_xz = sp.re(sp.diff(fun_phi, x, z))
        self.f_phi_yy = sp.re(sp.diff(fun_phi, y, y))
        self.f_phi_yz = sp.re(sp.diff(fun_phi, y, z))
        self.f_phi_zz = sp.re(sp.diff(fun_phi, z, z))

        self.f_elev = sp.re(fun_elev)
        self.f_elev_t = sp.re(sp.diff(fun_elev, t))

        self.f_elev_x = sp.re(sp.diff(fun_elev, x))
        self.f_elev_y = sp.re(sp.diff(fun_elev, y))

        self.f_elev_xx = sp.re(sp.diff(fun_elev, x, x))
        self.f_elev_xy = sp.re(sp.diff(fun_elev, x, y))
        self.f_elev_yy = sp.re(sp.diff(fun_elev, y, y))

    def Xfun(self, j):
        return sp.exp(- sp.I * j * self.dk * self.xswd)

    def Zfun_exact(self, j):
        kj = j * self.dk
        return sp.cosh(kj * (self.zswd + self.d)) / sp.cosh(kj * self.d)

    def Zhfun_exact(self, j):
        kj = j * self.dk
        return sp.sinh(kj * (self.zswd + self.d)) / sp.cosh(kj * self.d)

    def Zfun_exact_taylor(self, j):
        Z_e = self.Zfun_exact(j)
        Z_t = Z_e.series(z, 0, self.order).removeO()
        return sp.Piecewise((Z_e, self.zswd <= 0),
                            (Z_t, self.zswd > 0))

    def Zhfun_exact_taylor(self, j):
        Z_e = self.Zhfun_exact(j)
        Z_t = Z_e.series(z, 0, self.order).removeO()
        return sp.Piecewise((Z_e, self.zswd <= 0),
                            (Z_t, self.zswd > 0))

    def phi(self, x_app, y_app, z_app, t_app):
        return sp.N(self.f_phi.subs({x:x_app, y:y_app, z:z_app, t:t_app}))

    def stream(self, x_app, y_app, z_app, t_app):
        return sp.N(self.f_stream.subs({x: x_app, y: y_app, z: z_app, t: t_app}))

    def phi_t(self, x_app, y_app, z_app, t_app):
        return sp.N(self.f_phi_t.subs({x:x_app, y:y_app, z:z_app, t:t_app}))

    def grad_phi(self, x_app, y_app, z_app, t_app):
        phi_x = sp.N(self.f_phi_x.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_y = sp.N(self.f_phi_y.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_z = sp.N(self.f_phi_z.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        return {'x':phi_x, 'y':phi_y, 'z':phi_z}

    def grad_phi_2nd(self, x_app, y_app, z_app, t_app):
        phi_xx = sp.N(self.f_phi_xx.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_xy = sp.N(self.f_phi_xy.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_xz = sp.N(self.f_phi_xz.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_yy = sp.N(self.f_phi_yy.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_yz = sp.N(self.f_phi_yz.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_zz = sp.N(self.f_phi_zz.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        return {'xx':phi_xx, 'xy':phi_xy, 'xz':phi_xz,
                'yy':phi_yy, 'yz':phi_yz, 'zz':phi_zz}

    def acc_euler(self, x_app, y_app, z_app, t_app):
        phi_tx = sp.N(self.f_phi_tx.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_ty = sp.N(self.f_phi_ty.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
        phi_tz = sp.N(self.f_phi_tz.subs({x:x_app, y:y_app, z:z_app, t:t_app}))
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

    def dump_spectral_fun(self, j, dt, tmax):
        """For debug output plot file of time series for spectral component j"""
        out = open(os.path.join(self.tmpdir, 'spec_j%.dat' % j), 'w')
        out.write('t, h_j_re, h_j_im, ht_j_re, ht_j_im, c_j_re, c_j_im, ct_j_re, ct_j_im\n')
        t = 0.0
        T = sp.Symbol("T", real=True)  # SWD time
        hj = self.hfuns[j].fun(T)
        htj = sp.diff(hj, T)
        cj = self.cfuns[j].fun(T)
        ctj = sp.diff(cj, T)
        while t <= tmax:
            res_hj = hj.evalf(subs={T: t})
            res_htj = htj.evalf(subs={T: t})
            res_cj = cj.evalf(subs={T: t})
            res_ctj = ctj.evalf(subs={T: t})
            out.write('%f   %f %f   %f %f    %f %f   %f %f\n' % (t,
                       sp.re(res_hj), sp.im(res_hj),
                       sp.re(res_htj), sp.im(res_htj),
                       sp.re(res_cj), sp.im(res_cj),
                       sp.re(res_ctj), sp.im(res_ctj)))
            t += dt
        out.close()


    def write_swd(self, file_swd, dt, nsteps, too_short_file=False):

        out = open(file_swd, 'wb')
        out.write(pack('<f', 37.0221)) # Magic number
        out.write(pack('<i', 100))  # fmt
        out.write(pack('<i', 2))    # shp
        out.write(pack('<i', 1))    # amp
        out.write(pack('<30s', 'my_swd_symbolic'.ljust(30).encode('utf-8')))  # prog name
        out.write(pack('<20s', 'yyyy:mm:dd hh:ss'.ljust(20).encode('utf-8')))  # date
        wave_generator_data = "{'ole':1, 'dole':2, 'doffen':3}"
        nid = len(wave_generator_data)
        out.write(pack('<i', nid)) # length of input file
        out.write(pack('<{0}s'.format(nid), wave_generator_data.encode('utf-8')))   # Input file
        out.write(pack('<f', 9.81))  # acc. of gravity
        out.write(pack('<f', 1.0))   # lscale
        out.write(pack('<i', 0))     # nstrip
        out.write(pack('<i', nsteps))
        out.write(pack('<f', dt))
        out.write(pack('<i', self.order))
        out.write(pack('<i', self.n))
        out.write(pack('<f', self.dk))
        out.write(pack('<f', self.d))

        cf = []
        ctf = []
        hf = []
        htf = []
        T = sp.Symbol("T", real=True)  # SWD time
        for j in range(self.n + 1):
            hj = self.hfuns[j].fun(T)
            htj = sp.diff(hj, T)
            hf.append(hj)
            htf.append(htj)
            cj = self.cfuns[j].fun(T)
            ctj = sp.diff(cj, T)
            cf.append(cj)
            ctf.append(ctj)

        def dump(f, time):
            for j in range(self.n + 1):
                res = f[j].evalf(subs={T:time})
                out.write(pack('<f', sp.re(res)))
                out.write(pack('<f', sp.im(res)))

        if too_short_file:
            nout = nsteps // 2
        else:
            nout = nsteps
        for i in range(nout):
            t_swd = i*dt
            dump(hf, t_swd)
            dump(htf, t_swd)
            dump(cf, t_swd)
            dump(ctf, t_swd)

        out.close()

    def check_swd_meta(self, file_swd, n):
        if sys.version_info >= (3, 5):
            result = subprocess.run(["swd_meta", file_swd], capture_output=True)
            assert result.returncode == 0
            text = str(result.stdout)

            def check(tag, val):
                text_ok = "%-8s %s" % (tag + ':', val)
                assert text_ok in text, ("missing: %s" % text_ok)

            check("shp", 2)
            check("n", n)
