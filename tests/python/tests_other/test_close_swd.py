import os

import pytest

from spectral_wave_data import SpectralWaveData
from spectral_wave_data.tools import airy


@pytest.fixture(scope="module")
def make_swd(tmpdir_factory):
    mydir = str(tmpdir_factory.mktemp("close_tests"))
    file_swd = os.path.join(mydir, "my.swd")
    depth = 17.3
    grav = 9.81
    amps_inp = [1.2, 0.4, 0.8]
    dirs_inp = [173.2, 25.0, -130.0]
    phases_inp = [0.0, 210.0, 70.0]
    Twaves_inp = [3.0, 11.0, 70.0]
    airy.write_swd(file_swd, amps_inp, dirs_inp, phases_inp, Twaves=Twaves_inp, depth=depth, grav=grav)
    return file_swd


def test_close_close(make_swd):
    file_swd = make_swd
    swd = SpectralWaveData(file_swd)
    t = 3.0
    x = 1.4
    y = 4.5
    swd.update_time(t)
    zeta = swd.elev(x, y)
    swd.close()
    swd.close()


def test_close_update_time(make_swd):
    file_swd = make_swd
    swd = SpectralWaveData(file_swd)
    t = 3.0
    x = 1.4
    y = 4.5
    swd.update_time(t)
    zeta = swd.elev(x, y)
    swd.close()
    t *= 2
    with pytest.raises(AttributeError):
        swd.update_time(t)


def test_close_elev(make_swd):
    file_swd = make_swd
    swd = SpectralWaveData(file_swd)
    t = 3.0
    x = 1.4
    y = 4.5
    swd.update_time(t)
    zeta = swd.elev(x, y)
    swd.close()
    with pytest.raises(AttributeError):
        zeta_2 = swd.elev(x, y)


def test_garbage_collection(make_swd):
    file_swd = make_swd
    swd = SpectralWaveData(file_swd)
    t = 3.0
    x = 1.4
    y = 4.5
    swd.update_time(t)
    zeta = swd.elev(x, y)
    lst = [swd]
    del swd  # Only reduce reference count because object should still be accessible through lst[0]
    swd2 = lst[0]
    zeta2 = swd2.elev(x, y)
    assert zeta2 == pytest.approx(zeta)
    del lst  # Again, only reduce the reference count. swd2 is still alive
    zeta3 = swd2.elev(x, y)
    assert zeta3 == pytest.approx(zeta)
    swd2.close()
    with pytest.raises(AttributeError):
        zeta4 = swd2.elev(x, y)

