import numpy as np
import pytest
from spectral_wave_data import SpectralWaveData
import os


# make a list with the combination of the parameters to test
case_params = []
test_ids = []

shp_imp = [(1, 1), (2, 1), (4, 1), (4, 2), (5, 1)]
funcs = ['elev', 'grad_phi']
x0_y0_beta = [(0.0, 0.0, 0.0), (-23.67, 56.4, -23.65)]


z_ind = ['elev'] # functions independent of z

for shape, impl in shp_imp:
    for x0, y0, beta in x0_y0_beta:
        for nx in [32, 33, 34, 35]:  # hosm nx
            for ny in [1, 32, 33, 34, 35]: # hosm ny
                swd_name =  f'inputfiles/shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd'
                if not os.path.isfile(swd_name): # swd doesn't exist ==> illegal combination of shp, impl, nx, ny
                    continue
                for nsumx in [-1, 13, 14]:
                    for nsumy in [-1, 13, 14]:
                        for dc_bias in [True, False]:
                            for norder in [-1, 0, 2, 4, 7]:
                                for nx_fft in [-1, -2, 50, 51]:
                                    for ny_fft in [-1, -2, 50, 51]:
                                        # skip values for nsumy and ny_fft not valid for 1D
                                        if shape in [1, 2] and (nsumy > 1 or abs(ny_fft) != 1):
                                            continue
                                        for func in funcs:
                                            for z in [0, -2.1, 2.1]:
                                                # skip z-dependent for fields independent of z
                                                if func in z_ind and (norder != 0 or z != 0):
                                                    continue
                                                # skip norder for z <= 0
                                                if z <= 0 and norder != 0:
                                                    continue

                                                test_ids.append(f'{shape=}, {impl=}, {x0=}, {y0=}, {beta=}, {nx=}, {ny=}, {nsumx=}, {nsumy=}, {dc_bias=}, {norder=}, {nx_fft=}, {ny_fft=}, {func=}, {z=}')
                                                case_params.append((shape, impl, x0, y0, beta, nx, ny, nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft, func, z))


@pytest.mark.parametrize('test_case', case_params, ids=test_ids)  # tests to run through
def test_fft(test_case):
    shape, impl, x0, y0,beta, nx, ny, nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft, func, z = test_case

    with SpectralWaveData(f'inputfiles/shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd',
                          x0=x0, y0=y0, beta=beta,
                          nsumx=nsumx, nsumy=nsumy, norder=norder, dc_bias=dc_bias) as swd:
        swd.update_time(0.0)

        x_fft, y_fft = swd.xy_fft(nx_fft=nx_fft, ny_fft=ny_fft)
        assert(x_fft.shape == y_fft.shape)

        # test only nxp x nyp gridpoints
        nxp, nyp = 5, 5
        ixs = np.unique(np.linspace(0, x_fft.shape[0] - 1, nxp).astype(int))
        iys = np.unique(np.linspace(0, y_fft.shape[1] - 1, nyp).astype(int))

        if func in z_ind:
            arr_fft = getattr(swd, f'{func}_fft')(nx_fft=nx_fft, ny_fft=ny_fft)
        else:
            arr_fft = getattr(swd, f'{func}_fft')(z, nx_fft=nx_fft, ny_fft=ny_fft)

        # scalar output
        if arr_fft.ndim == 2:
            res = np.zeros((ixs.size, iys.size))
            res_fft = np.zeros((ixs.size, iys.size))
            assert(x_fft.shape == arr_fft.shape)
            for i, ix in enumerate(ixs):
                for j, iy in enumerate(iys):
                    if func in z_ind:
                        res[i, j] = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy])
                    else:
                        res[i, j] = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy], z)
                    res_fft[i, j] = arr_fft[ix, iy]
            assert(np.allclose(res, res_fft))
        elif arr_fft.ndim == 3: # 3-component output
            res = np.zeros((3, ixs.size, iys.size))
            res_fft = np.zeros((3, ixs.size, iys.size))
            assert((3, x_fft.shape[0], x_fft.shape[1]) == arr_fft.shape)
            for i, ix in enumerate(ixs):
                for j, iy in enumerate(iys):
                    if func in z_ind:
                        resxyz = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy])
                    else:
                        resxyz = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy], z)
                    res[:, i, j] = np.asarray([resxyz.x, resxyz.y, resxyz.z])
                    res_fft[:, i, j] = arr_fft[:, ix, iy]
            print(np.max(np.abs(res[0, :, :] - res_fft[0, :, :])))
            print(np.max(np.abs(res[1, :, :] - res_fft[1, :, :])))
            print(np.max(np.abs(res[2, :, :] - res_fft[2, :, :])))
            assert(np.allclose(res[0, :, :], res_fft[0, :, :]))
            assert(np.allclose(res[1, :, :], res_fft[1, :, :]))
            assert(np.allclose(res[2, :, :], res_fft[2, :, :]))


