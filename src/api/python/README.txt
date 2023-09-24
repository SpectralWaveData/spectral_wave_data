This directory is for building the spectral_wave_data library (whl install file) for Python.

You need a recent Fortran and C compiler (later than 2017) to build this.
Precompiled libraries can be downloaded from the Github repository 'spectral_wave_data_precompiled'

---------------------------------------------------------------------------------------------
Do the following steps on Windows:
---------------------------------------------------------------------------------------------
1) cd Cmake
2) Adjust the actual CMakeBuild*.bat file to reflect your actual version of Visual Studio.
3) Run this CMakeBuild*.bat to create the the relevant Visual Studio Solution (32 or 64bit)
4) cd actual Build_* folder.
5) Open spectral_wave_data.sln with Visual Studio
6) Make sure the 'Release' flag and actual Binary configuration (64 vs 32bit) is selected in the tool bar of Visual Studio.
7) Build -> Build Solution
8) Check that the solution build without error messages and check that SpectralWaveData.dll is created in the Release sub-folder.
9) Enter the Python/distribution folder.
10) Edit setup.py to set correct version number for the new distribution.
11) Run the make_wheel_*.bat from a terminal using the Python version you want to create support for.
    You need to be in a 64-bit Python environment to run the *_win64.bat version
    You need to be in a 32-bit Python environment to run the *_win32.bat version
12) If no errors a *.whl file is copied to the python/whl folder ready for distribution.
13) Test new whl file: 
       pip uninstall spectral_wave_data (if already installed)
       pip install spectral_wave_data-*.whl
       cd spectral-wave-data/src/api_examples
       run application_swd_api.py
14) Run cleanup.bat (Also between making 32 and 64 bit versions)

---------------------------------------------------------------------------------------------
Do the following steps on Linux:
---------------------------------------------------------------------------------------------
1) run BuildWheelsLinux.bash to build wheels for both python 2.7 and python 3.
2) Test new whl file(s): 
       pip uninstall spectral_wave_data (if already installed)
       pip install spectral_wave_data-*.whl
       cd spectral-wave-data/src/api_examples
       run application_swd_api.py

If BuildWheelsLinux.bash fails, it is recommended to carry out the procedure in the following steps to identify the problem:
1) cd Cmake
2) Run relevant version of CMakeBuild_Py*_Linux.bash
3) Make sure libSpectralWaveData.so is created in the Build_Py*_Linux folder
4) Enter the Python/distribution folder.
5) If needed: Edit setup.py to set correct version number for the new distribution. (E.g. 2.0.6)
6) Possibly edit make_wheel_*_linux_*.bash path to the correct python binary in line 11.
7) Run make_wheel_*_linux_*.bash
8) If no errors a *.whl file is copied to the Python/whl folder ready for distribution.
9) Test new whl file: 
       pip uninstall hosm (if already installed)
       pip install spectral_wave_data-*.whl
       cd spectral-wave-data/src/api_examples
       run application_swd_api.py
10) Run cleanup.bash


Jens B. Helmers, Odin Gramstad
DNVGL, 2019
