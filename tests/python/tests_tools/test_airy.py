import pytest
import numpy as np

from spectral_wave_data.tools import airy

def test_omega2kwave():
    depth = 10.0
    grav = 9.81
    Cvals =[0.0, 0.0001, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.3, 2.5, 2.8,
            3.0, 3.2, 3.5, 4.0, 10.0, 20.0, 50000.0]
    for c in Cvals:
        omega = np.sqrt(c * grav / depth)
        kwave = airy.omega2kwave(omega, depth, grav)
        omega2 = kwave * grav * np.tanh(kwave * depth)
        assert abs(omega ** 2 - omega2) <= 100 * np.finfo(omega).eps
    # infinte water
    depth = -1.0
    omega = 2.71
    kwave = airy.omega2kwave(omega, depth, grav)
    omega2 = kwave * grav
    assert abs(omega ** 2 - omega2) <= 100 * np.finfo(omega).eps



