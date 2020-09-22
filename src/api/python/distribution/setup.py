# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import platform
import os
import glob
import sys

from distutils.util import get_platform
from setuptools import setup, find_packages, Distribution
from setuptools.command.install import install
from setuptools.command.test import test as TestCommand
from wheel.bdist_wheel import bdist_wheel


# To read version number from version_definition.py in root directory
sys.path.insert(0, os.path.abspath(os.path.join("..", "..", "..", "..")))
from version_definition import version_full

# Get the long description from 'readme_pypi.md'
# with open(os.path.join('readme_pypi.md'), encoding='utf-8') as f:  # Not legal in py2.7
with open(os.path.join("readme_pypi.md")) as f:
    long_description = f.read()


# Copy files into the package directory
pkg_dir = "spectral_wave_data"
if platform.system() == "Windows":
    libfiles = [
        os.path.join(pkg_dir, "SpectralWaveData.dll"),
        os.path.join(pkg_dir, "SpectralWaveData.exp"),
        os.path.join(pkg_dir, "spectral_wave_data.h"),
        os.path.join(pkg_dir, "SpectralWaveData.lib"),
        os.path.join(pkg_dir, "SWD_logo.ico"),
        os.path.join(pkg_dir, "README.txt"),
        os.path.join(pkg_dir, "LICENSE.txt"),
    ]
    libfiles += glob.glob(os.path.join(pkg_dir, "*.mod"))
elif platform.system() == "Linux":
    libfiles = [
        os.path.join(pkg_dir, "libSpectralWaveData.so"),
        os.path.join(pkg_dir, "spectral_wave_data.h"),
        os.path.join(pkg_dir, "SWD_logo.ico"),
        os.path.join(pkg_dir, "README.txt"),
        os.path.join(pkg_dir, "LICENSE.txt"),
    ]
    # Not required with ifort and intel static compilation
    libfiles += glob.glob(os.path.join(pkg_dir, "*.mod"))
else:
    raise AssertionError("Not supported platform: " + platform.system())

out = open("MANIFEST.in", "w")
for f in libfiles:
    out.write("include {0}\n".format(f))
out.close()


########################################################################
# Custom commands

# Make the setup.py test command work
class PyTest(TestCommand):
    description = "Run spectral_wave_data's tests with pytest"

    def initialize_options(self):
        TestCommand.initialize_options(self)

    def finalize_options(self):
        TestCommand.finalize_options(self)

    def run_tests(self):
        import pytest

        args = ["-v", "--durations=10"]
        if self.verbose:
            args.append("-s")
        here = os.path.dirname(__file__)
        args.append(os.path.join(here, "tests/python/"))
        errno = pytest.main(args)
        sys.exit(errno)


class my_install(install):
    """
    Force use of "platlib" dir for auditwheel to recognize this
    as a non-pure build, needed for manylinux support
    """

    def finalize_options(self):
        install.finalize_options(self)
        self.install_libbase = self.install_platlib
        self.install_lib = self.install_platlib


class my_bdist_wheel(bdist_wheel):
    """
    Ensure we are marked as a universal, binary distribution.
    Not dependent on Python version, but dependent on platform
    """

    def finalize_options(self):
        bdist_wheel.finalize_options(self)
        # Mark us as not a pure python package
        self.root_is_pure = False
        self.universal = True
        self.plat_name_supplied = True
        self.plat_name = get_platform()

    def get_tag(self):
        python, abi, plat = bdist_wheel.get_tag(self)
        # The Python source (and the compiled library for ctypes)
        # does work on both Py2 and Py3, so no need to constrain
        # anthing but the platform (due to the compiled library)
        python, abi = "py2.py3", "none"
        return python, abi, plat


class BinaryDistribution(Distribution):
    """
    Distribution which always forces a binary package with platform name
    Uses the my_install and my_bdist_wheel commands from above
    """

    def __init__(self, *args, **kwargs):
        Distribution.__init__(self, *args, **kwargs)
        self.cmdclass["install"] = my_install
        self.cmdclass["bdist_wheel"] = my_bdist_wheel

    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False


########################################################################
# The Python package configuration

setup(
    name="spectral_wave_data",
    version=version_full,
    description="An API for ocean wave kinematics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # The project's main homepage.
    url="https://github.com/SpectralWaveData",
    maintainer="Jens B. Helmers & Odin Gramstad, DNVGL",
    maintainer_email="Jens.Bloch.Helmers@dnvgl.com",
    license="MIT license - Copyright DNVGL 2019-2020",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: MIT License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        # Supported operating systems
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
    ],
    # What does your project relate to?
    keywords="API ocean wave kinematics",
    packages=find_packages(),
    setup_requires=["wheel"],
    entry_points={"console_scripts": ["swd_meta=spectral_wave_data.scripts.swd_meta:main"],},
    include_package_data=True,
    # install_requires=['numpy'],
    # Current testing is too time consuming for automatic testing!!!
    # Configure the "test" command
    # tests_require=['numpy', 'sympy', 'pytest', 'pytest-cov', 'pytest-html', 'pytest-codestyle'],
    # cmdclass={'test': PyTest},
    distclass=BinaryDistribution,
)
