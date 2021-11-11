import os
import numpy as np
import pytest
from spectral_wave_data import SpectralWaveData

input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "inputfiles")

# make a list with the combination of the parameters to test
case_params = []
test_ids = []
shp_imp = [(1, 1), (2, 1), (4, 1), (4, 2), (5, 1)]
funcs = ['elev', 'elev_t', 'phi', 'phi_t']
for shape, impl in shp_imp:
    for func in funcs:
      test_ids.append(f'{shape=}, {impl=}, {func=}')
      case_params.append((shape, impl, func))

@pytest.mark.parametrize('test_case', case_params, ids=test_ids)  # tests to run through
def test_onetimestep(test_case):
    shape, impl, func = test_case

    nstep_max = 5 # swd file with the largest number of steps
    dt = 10.0  # delta t in the swd files

    for ut in range(nstep_max): # check time step ut in all .swd files with different nsteps
        res = []
        for nsteps in range(ut+1, nstep_max+1):
            file_swd = os.path.join(input_dir, f"shape{shape}_impl{impl}_nsteps{nsteps}.swd")
            with SpectralWaveData(file_swd) as swd:
                for i in range(10): # randomly jump a bit chack and forth in time
                    t_rand = np.random.uniform(0.0, (nsteps - 1)*dt)
                    swd.update_time(t_rand)
                    # update for time ut
                    swd.update_time(ut*dt)
                    if 'elev' in func:
                        res.append(getattr(swd, func)(1.23, 4.32))
                    else:
                        res.append(getattr(swd, func)(1.23, 4.32, 0.34))
        assert(np.allclose(res, res[0])) # make sure all .swd files give the same result
