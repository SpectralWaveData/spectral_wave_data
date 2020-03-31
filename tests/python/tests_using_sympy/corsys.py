from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import sympy as sp

from symbols import x, y, z, t

class CorSys(object):

    def __init__(self, x0, y0, t0, beta):
        self.x0 = x0
        self.y0 = y0
        self.t0 = t0
        self.beta = beta
        self.cbeta = sp.cos(beta * sp.pi / 180)
        self.sbeta = sp.sin(beta * sp.pi / 180)

    def app2swd(self):
        xswd = self.x0 + x * self.cbeta + y * self.sbeta
        yswd = self.y0 - x * self.sbeta + y * self.cbeta
        zswd = z
        tswd = self.t0 + t
        return xswd, yswd, zswd, tswd

