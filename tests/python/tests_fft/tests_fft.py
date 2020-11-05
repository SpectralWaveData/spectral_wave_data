import numpy as np
import pytest
from spectral_wave_data import SpectralWaveData
import os


# make a list with the combination of the parameters to test
case_params = []
test_ids = []

# shp_imp = [(1, 1), (2, 1), (4, 1), (4, 2), (5, 1)]
shp_imp = [(1, 1), (2, 1), (4, 1), (5, 1)]


z_ind = ['elev'] # functions independent of z

for shape, impl in shp_imp:
    for nx in [32, 33, 34, 35]:  # hosm nx
        for ny in [1, 32, 33, 34, 35]: # hosm ny
            swd_name =  f'inputfiles/shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd'
            if not os.path.isfile(swd_name): # swd doesn't exist ==> illegal combination of shp, impl, nx, ny
                continue
            for nsumx in [-1, 13, 14]:
                for nsumy in [-1, 13, 14]:
                    for dc_bias in [True, False]:
                        for norder in [0, -1]:
                            for nx_fft in [-1, -2, 50, 51]:
                                for ny_fft in [-1, -2, 50, 51]:
                                    # skip values for nsumy and ny_fft not valid for 1D
                                    if shape in [1, 2] and (nsumy > 1 or abs(ny_fft) != 1):
                                        continue
                                    for func in ['elev']:
                                        for z in [0, -1.1, 1.1]:
                                            # skip z-dependent for fields independent of z
                                            if func in z_ind and (norder != 0 or z != 0):
                                                continue

                                            test_ids.append(f'{shape=}, {impl=}, {nx=}, {ny=}, {nsumx=}, {nsumy=}, {dc_bias=}, {norder=}, {nx_fft=}, {ny_fft=}, {func=}, {z=}')
                                            case_params.append((shape, impl, nx, ny, nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft, func, z))


@pytest.mark.parametrize('test_case', case_params, ids=test_ids)  # tests to run through
def test_fft(test_case):
    shape, impl, nx, ny, nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft, func, z = test_case

    swd = SpectralWaveData(f'inputfiles/shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd', nsumx=nsumx, nsumy=nsumy,
                           norder=norder, dc_bias=dc_bias)
    swd.update_time(0.0)

    x_fft = swd.x_fft(nx_fft=nx_fft)
    y_fft = swd.y_fft(ny_fft=ny_fft)

    # test only nxp x nyp gridpoints
    nxp, nyp = 5, 5
    ixs = np.unique(np.linspace(0, x_fft.size - 1, nxp).astype(int))
    iys = np.unique(np.linspace(0, y_fft.size - 1, nyp).astype(int))
    res = np.zeros((ixs.size, iys.size))
    res_fft = np.zeros((ixs.size, iys.size))

    arr_fft = getattr(swd, f'{func}_fft')(nx_fft=nx_fft, ny_fft=ny_fft)
    assert((x_fft.size, y_fft.size) == arr_fft.shape)
    for i, ix in enumerate(ixs):
        for j, iy in enumerate(iys):
            res[i, j] = getattr(swd, func)(x_fft[ix], y_fft[iy])
            res_fft[i, j] = arr_fft[ix, iy]
    swd.close()
    assert(np.allclose(res, res_fft))

