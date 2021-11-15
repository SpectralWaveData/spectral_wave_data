import os
import sys

import pytest

from spectral_wave_data import SpectralWaveData, SwdError
from spectral_wave_data.tools import airy




@pytest.fixture(scope="module")
def make_swd(tmpdir_factory):
    mydir = str(tmpdir_factory.mktemp("with_test"))
    file_swd = os.path.join(mydir, "my.swd")
    depth = 17.3
    grav = 9.81
    amps_inp = [1.2, 0.4, 0.8]
    dirs_inp = [173.2, 25.0, -130.0]
    phases_inp = [0.0, 210.0, 70.0]
    Twaves_inp = [3.0, 11.0, 70.0]
    airy.write_swd(file_swd, amps_inp, dirs_inp, phases_inp, Twaves=Twaves_inp, depth=depth, grav=grav)
    return file_swd


def test_with_1(make_swd):
    file_swd = make_swd
    with SpectralWaveData(file_swd) as swd:
        assert swd["n"] == 3


def test_with_2(make_swd):
    file_swd = make_swd
    beta = 53.2
    t = 5.4
    x = 3.3
    y = 4.4
    # Do it using the with statement...
    with SpectralWaveData(file_swd, beta=beta) as swd:
        swd.update_time(t)
        zeta_1 = swd.elev(x, y)
    # Do it the classical way...
    swd = SpectralWaveData(file_swd, beta=beta)
    swd.update_time(t)
    zeta_2 = swd.elev(x, y)
    swd.close()
    assert zeta_1 == pytest.approx(zeta_2)


@pytest.mark.skipif(sys.platform=="linux", reason="pytest on Linux has problem raising SwdError, not python")
def test_with_3(make_swd):
    file_swd = make_swd
    with SpectralWaveData(file_swd) as swd:
        beta = swd["beta"]
        with pytest.raises(SwdError):
            x = swd["asfasdf"]


# This test works! However:
# pytest has some general issues with OSError. Just ignore the clutter screen dump in output....
#def test_with_4(make_swd):
#    file_swd = make_swd
#    with SpectralWaveData(file_swd) as swd:
#        beta = swd["beta"]
#
#    # SWD object should be closed by now...
#    with pytest.raises(OSError):
#        x0 = swd["x0"]
