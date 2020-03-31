# -*- coding: utf-8 -*-
   
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

from setuptools import setup, Distribution, find_packages
from setuptools.command.test import test as TestCommand
import platform, os, glob, sys

# To read version number from version_definition.py in root directory
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', '..', '..')))
from version_definition import version_full

# Get the long description from 'readme_pypi.md'
#with open(os.path.join('readme_pypi.md'), encoding='utf-8') as f:  # Not legal in py2.7
with open(os.path.join('readme_pypi.md')) as f:
    long_description = f.read()

# Make the setup.py test command work
class PyTest(TestCommand):
    description = 'Run spectral_wave_data\'s tests with pytest'

    def initialize_options(self):
        TestCommand.initialize_options(self)

    def finalize_options(self):
        TestCommand.finalize_options(self)

    def run_tests(self):
        import pytest
        args = ['-v', '--durations=10']
        if self.verbose:
            args.append('-s')
        args.append(os.path.join(here, 'tests/python/'))
        errno = pytest.main(args)
        sys.exit(errno)


if platform.system()=='Windows':
    dir = 'spectral_wave_data'
    libfiles = [os.path.join(dir, 'SpectralWaveData.dll'),
                os.path.join(dir, 'SpectralWaveData.exp'),
                os.path.join(dir, 'spectral_wave_data.h'),
                os.path.join(dir, 'SpectralWaveData.lib'),
                os.path.join(dir, 'SWD_logo.ico'),
                os.path.join(dir, 'README.txt'),
                os.path.join(dir, 'LICENSE.txt')]
    libfiles += glob.glob(os.path.join(dir, '*.mod'))
elif platform.system()=='Linux':
    dir = 'spectral_wave_data'
    libfiles = [os.path.join(dir, 'libSpectralWaveData.so'),
                os.path.join(dir, 'spectral_wave_data.h'),
                os.path.join(dir, 'SWD_logo.ico'),
                os.path.join(dir, 'README.txt'),
                os.path.join(dir, 'LICENSE.txt')]
    libfiles += glob.glob(os.path.join(dir, '*.mod'))
else:
    raise AssertionError('Not supported platform: ' + platform.system())

out = open('MANIFEST.in', 'w')
for f in libfiles:
    out.write('include {0}\n'.format(f))
out.close()

class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name"""
    def has_ext_modules(foo):
        return True

setup(
    name='spectral_wave_data',
    version=version_full,

    description='An API for ocean wave kinematics',
    long_description=long_description,
    long_description_content_type='text/markdown',

    # The project's main homepage.
    url='https://github.com/SpectralWaveData',

    maintainer='Jens B. Helmers & Odin Gramstad, DNVGL',
    maintainer_email='Jens.Bloch.Helmers@dnvgl.com',

    license='MIT license - Copyright DNVGL 2019-2020',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',

        # Supported operating systems
        'Operating System :: Microsoft :: Windows :: Windows 10',
        'Operating System :: POSIX :: Linux'

    ],

    # What does your project relate to?
    keywords='API ocean wave kinematics',

    packages=find_packages(),
    setup_requires=['wheel'],
    entry_points = {
        'console_scripts': ['swd_meta=spectral_wave_data.scripts.swd_meta:main'],
    },
    include_package_data=True,
    #install_requires=['numpy'],

    # Current testing is too time consuming for automatic testing!!!
    # Configure the "test" command
    #tests_require=['numpy', 'sympy', 'pytest', 'pytest-cov', 'pytest-html', 'pytest-codestyle'],
    #cmdclass={'test': PyTest},

    distclass=BinaryDistribution
)
