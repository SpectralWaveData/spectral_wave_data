import os
import sys
import subprocess
from struct import pack

import pytest

from spectral_wave_data import SpectralWaveData, SwdAllocateError, SwdFileCantOpenError, \
    SwdFileDataError, SwdFileBinaryError, SwdInputValueError
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


def test_open_not_exist(tmp_path):
    file_swd = os.path.join(str(tmp_path), "this_file_does_not_exist.swd")
    with pytest.raises(SwdFileCantOpenError):
        swd = SpectralWaveData(file_swd)


def test_open_empty(tmp_path):
    file_swd = os.path.join(str(tmp_path), "empty.swd")
    out = open(file_swd, 'wb')
    out.close()
    with pytest.raises(SwdFileDataError):
        swd = SpectralWaveData(file_swd)


def test_open_double_precision(tmp_path):
    file_swd = os.path.join(str(tmp_path), "my_double_precision.swd")
    out = open(file_swd, 'wb')
    out.write(pack('<d', 37.0221))  # Magic number
    out.close()
    with pytest.raises(SwdFileBinaryError):
        swd = SpectralWaveData(file_swd)


def test_open_big_endian(tmp_path):
    file_swd = os.path.join(str(tmp_path), "my_big_endian.swd")
    out = open(file_swd, 'wb')
    out.write(pack('>f', 37.0221))  # Magic number
    out.close()
    with pytest.raises(SwdFileBinaryError):
        swd = SpectralWaveData(file_swd)


def test_open_big_endian_double_precision(tmp_path):
    file_swd = os.path.join(str(tmp_path), "my_big_endian_double_precision.swd")
    out = open(file_swd, 'wb')
    out.write(pack('>d', 37.0221))  # Magic number
    out.close()
    with pytest.raises(SwdFileBinaryError):
        swd = SpectralWaveData(file_swd)


def test_open_not_a_swd_file(tmp_path):
    file_swd = os.path.join(str(tmp_path), "my_not_swd.swd")
    out = open(file_swd, 'w')
    out.write("This is not a SWD file")
    out.close()
    with pytest.raises(SwdFileDataError):
        swd = SpectralWaveData(file_swd)


def test_data_in_swd(tmp_path):
    file_swd = os.path.join(str(tmp_path), "bad_fmt_value.swd")
    out = open(file_swd, 'wb')
    out.write(pack('<f', 37.0221))  # Magic number
    out.write(pack('<i', 100000))   # fmt value
    out.close()
    with pytest.raises(SwdFileDataError):
        swd = SpectralWaveData(file_swd)   # fmt value in file is too high


def test_bad_constructor_input_value(make_swd):
    file_swd = make_swd
    with pytest.raises(SwdInputValueError):
        swd = SpectralWaveData(file_swd, norder=1000)  # Not allowed for Airy waves
