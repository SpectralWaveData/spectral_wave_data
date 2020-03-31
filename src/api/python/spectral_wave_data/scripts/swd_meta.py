from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import sys

from spectral_wave_data import SpectralWaveData, SwdFileCantOpenError, \
     SwdFileBinaryError, SwdFileDataError, SwdInputValueError, \
     SwdAllocateError

def main():

    narg = len(sys.argv)
    if narg != 2:
        print("Usage: swd_meta my.swd")
        sys.exit()
    file_swd = sys.argv[-1]

    try:
        swd = SpectralWaveData(file_swd, x0=0.0, y0=0.0, t0=0.0, beta=0.0)
    except SwdFileCantOpenError as e:
        print("Not able to open: %s" % file_swd)
        sys.exit()
    except SwdFileBinaryError as e:
        print("This SWD file don't have the correct binary convention: %s" % file_swd)
        sys.exit()
    except SwdFileBinaryError as e:
        print("This file don't look like a SWD-file: %s" % file_swd)
        sys.exit()

    def write(tag):
        print('%-8s %s' % (tag + ':', swd[tag]))

    write('version')
    write('prog')
    write('date')
    write('fmt')
    write('shp')
    write('amp')
    write('tmax')
    write('dt')
    write('nsteps')
    write('nstrip')
    write('order')
    write('d')

    shp = swd['shp']
    if shp in [1, 2, 3]:
        # Long-crested seas
        write('n')
        if shp == 3:
            write('nh')
        write('sizex')
        write('lmax')
        write('lmin')
        write('dk')

    if shp in [4, 5]:
        # Short-crested seas
        write('nx')
        write('ny')
        write('sizex')
        write('sizey')
        write('lmax')
        write('lmin')
        write('dkx')
        write('dky')

    if shp in [6]:
        # Airy waves
        write('n')

    write('cid')
