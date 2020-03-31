"""
This script converts a vector graphic SVG file to a multi-resolution Windows *.ico file and a set of png files.

You need to install:
 1) Inkscape (inkscape.org) for a very good SVG to PNG engine.
 2) ImageMagick (imagemagick.org) for converting multiple png into a multiple resolution *.ico file

"""
import os

svg = 'SWD_logo.svg'
svg_tiny = 'SWD_logo_tiny.svg'  # No letters are visible on 16x16 pixels

path_inckscape = r'C:\Program Files\Inkscape\inkscape.exe'
path_imagemagick = r'C:\Program Files\ImageMagick-7.0.9-Q16\magick.exe'


# Requested resolutions
res = [16, 32, 48, 64, 128, 256, 512, 1024, 2048]

pngs = ''
for r in res:
    png = '%s_%ix%i.png' % (svg[:-4], r, r)
    if r < 500:   # Very large png should not be included in the ico file.
        pngs += png + ' '
    if r > 16:
        cmd = r'"%s" -z -e %s -w %i -h %i %s' % (path_inckscape, png, r, r, svg)
    else:
        cmd = r'"%s" -z -e %s -w %i -h %i %s' % (path_inckscape, png, r, r, svg_tiny)
    os.system(cmd)

cmd = r'"%s" %s -colors 256 %s' % (path_imagemagick, pngs, svg[:-4] + '.ico')
os.system(cmd)
