from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import copy

import sympy as sp
import numpy as np


class Tfun(object):

    def __init__(self, poly_order):
        self.porder = poly_order
        self.c = []
        for i in range(self.porder + 1):
            # The swd file apply float which has less precision than python.
            # Consequently, we only pick a few digits for the random
            # temporal coefficients.
            re = np.around(np.random.uniform(-1.0, 1.0), decimals=4)
            im = np.around(np.random.uniform(-1.0, 1.0), decimals=4)
            self.c.append(re + sp.I * im)

    def fun(self, tswd):
        res = 0
        for i in range(self.porder + 1):
            res += self.c[i] * tswd**i
        return res

    def scale(self, fac):
        res = copy.deepcopy(self)
        for i in range(len(self.c)):
            res.c[i] = self.c[i] * fac
        return res

def create_Tfuns_1d(n_swd, ipol=0, seed=0):
    """
    Define all temporal functions using polynomial of proper order
    related to the the ipol-scheme applied in the numerical method.
    Consequently, all temporal calculations should be exact.
    The polynomial coefficients are random in [-1.0, 1.0]

    :param n_swd:   The number of non-zero shape functions in model
    :param ipol:   The Fortran/C/C++/Python constructor parameter
                    0 = Default (C^2 continous scheme)
    """

    if ipol == 0:
        poly_order = 5
    elif ipol == 1:
        poly_order = 3
    elif ipol == 2:
        poly_order = 7
    else:
        raise("Not expected value of ipol!!!")

    np.random.seed(seed)  # All polynomial coefficients are random in [-1.0, 1.0]
    cfuns = []
    hfuns = []
    for j in range(n_swd + 1):
        cfuns.append(Tfun(poly_order))
        hfuns.append(Tfun(poly_order))
    return cfuns, hfuns


def create_Tfuns_2d(nx_swd, ny_swd, ipol=0, seed=0):
    """
    Define all temporal functions using polynomial of proper order
    related to the the ipol-scheme applied in the numerical method.
    Consequently, all temporal calculations should be exact.
    The polynomial coefficients are random in [-1.0, 1.0]

    :param nx_swd:   The number of non-zero shape functions in x-dir

    :param ny_swd:   The number of non-zero shape functions in y-dir

    :param ipol:   The Fortran/C/C++/Python constructor parameter
                    0 = Default (C^2 continous scheme)
                    1 = C^1 continous
                    2 = C^3 continous

    :param seed:   Random number generator seed

    :return: lists of Tfun for c and h functions

    """

    if ipol == 0:
        poly_order = 5
    elif ipol == 1:
        poly_order = 3
    elif ipol == 2:
        poly_order = 7
    else:
        raise("Not expected value of ipol!!!")

    np.random.seed(seed)  # All polynomial coefficients are random in [-1.0, 1.0]
    cfuns = np.empty((nx_swd + 1, 2 * ny_swd + 1), dtype=Tfun)
    hfuns = np.empty((nx_swd + 1, 2 * ny_swd + 1), dtype=Tfun)
    for jx in range(nx_swd + 1):
        for jy in range(-ny_swd, ny_swd + 1):
            cfuns[jx, ny_swd + jy] = Tfun(poly_order)
            hfuns[jx, ny_swd + jy] = Tfun(poly_order)
    return cfuns, hfuns