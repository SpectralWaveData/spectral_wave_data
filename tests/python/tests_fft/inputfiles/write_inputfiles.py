#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:23:52 2020

@author: oding
"""
import f90nml
import subprocess
import os

def write_inputfile(shape, impl, nx, ny):

    out_name = f'shape{shape}_impl{impl}_nx{nx}_ny{ny}'
    inp_name = f'{out_name}.inp'

    nml = dict()
    nml['Hs'] = 7.6
    nml['Tp'] = 13.0
    if shape in [1, 4]:
        nml['h'] = -1
    elif shape in [2, 5]:
        nml['h'] = 35.0
    if shape in [1, 2]:
        if ny != 1: return
        nml['ndim'] = 1
        nml['nx'] = nx
    elif shape in [4, 5]:
        if ny == 1: return
        nml['ndim'] = 2
        nml['nx'] = nx
        nml['ny'] = ny
        if shape == 4 and impl == 1 and nx == ny: return
        if shape == 4 and impl == 2 and nx != ny: return

    nml['tmax'] = 0.0
    nml['sf_nout'] = 1
    nml['outfile'] = out_name
    nml['output'] = ['swd']

    f90nml.Namelist(hosm=nml).write(inp_name, force=True)
    return inp_name

shp_imp = [(1, 1), (2, 1), (4, 1), (4, 2), (5, 1)]

exe = '/home/oding/DNV/prog_dev/Wamod/HOSM/bin/HOSM_time_series_nxny'
for shape, impl in shp_imp:
    for nx in [32, 33, 34, 35]:
        for ny in [1, 32, 33, 34, 35]:
            inp_name = write_inputfile(shape, impl, nx, ny)
            if inp_name is None: continue
            subprocess.run(f'{exe} {inp_name}', shell=True, env={"LD_LIBRARY_PATH": "/usr/local/hdf5/lib"})
            os.remove(inp_name)
            os.remove(inp_name.replace('.inp', '.h5'))
