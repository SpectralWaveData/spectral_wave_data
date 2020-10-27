import os

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


def test_get_components_finite_depth(tmp_path):

    file_swd = os.path.join(str(tmp_path), "my_shp6.swd")

    depth = 17.3
    grav = 9.81
    nwaves = 3
    amps_inp = [1.2, 0.4, 0.8]
    dirs_inp = [173.2, 25.0, -130.0]
    phases_inp = [0.0, 210.0, 70.0]
    Twaves_inp = [3.0, 11.0, 70.0]

    airy.write_swd(file_swd, amps_inp, dirs_inp, phases_inp, Twaves=Twaves_inp, depth=depth, grav=grav)
    res = airy.get_components(file_swd)

    assert res["n"] == nwaves
    assert res["depth"]  == pytest.approx(depth)
    assert res["grav"] == pytest.approx(grav)
    for i in range(nwaves):
        assert res["amps"][i] == pytest.approx(amps_inp[i])
        assert res["dirs_deg"][i] == pytest.approx(dirs_inp[i])
        assert res["phases_deg"][i] == pytest.approx(phases_inp[i])
        assert res["Twaves"][i] == pytest.approx(Twaves_inp[i])
        omg_inp = 2.0 * np.pi / Twaves_inp[i]
        assert res["omegas"][i] == pytest.approx(omg_inp)
        kw_inp = airy.omega2kwave(omg_inp, depth, grav)
        assert res["kwaves"][i] == pytest.approx(kw_inp)
        lw_inp = 2.0 * np.pi / kw_inp
        assert res["Lwaves"][i] == pytest.approx(lw_inp)


def test_get_components_infinite_depth(tmp_path):

    file_swd = os.path.join(str(tmp_path), "my_shp6_inf.swd")

    depth = -1.0
    grav = 9.81
    nwaves = 1
    amps_inp = [10.2]
    dirs_inp = [17.2]
    phases_inp = [-10.0]
    Twaves_inp = [13.0]

    airy.write_swd(file_swd, amps_inp, dirs_inp, phases_inp, Twaves=Twaves_inp, depth=depth, grav=grav)
    res = airy.get_components(file_swd)

    assert res["n"] == nwaves
    assert res["depth"]  == pytest.approx(depth)
    assert res["grav"] == pytest.approx(grav)
    for i in range(nwaves):
        assert res["amps"][i] == pytest.approx(amps_inp[i])
        assert res["dirs_deg"][i] == pytest.approx(dirs_inp[i])
        assert res["phases_deg"][i] == pytest.approx(phases_inp[i])
        assert res["Twaves"][i] == pytest.approx(Twaves_inp[i])
        omg_inp = 2.0 * np.pi / Twaves_inp[i]
        assert res["omegas"][i] == pytest.approx(omg_inp)
        kw_inp = airy.omega2kwave(omg_inp, depth, grav)
        assert res["kwaves"][i] == pytest.approx(kw_inp)
        lw_inp = 2.0 * np.pi / kw_inp
        assert res["Lwaves"][i] == pytest.approx(lw_inp)
