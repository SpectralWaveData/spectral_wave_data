This directory is for building a C++ test application.
First it builds the static library spectral-wave-data which is linked into the program.

The kind_values.f90 file is important to ensure that the numerical precision
(float or double) applied in the C++-program match the fortran code.

If the C++-program apply 'double' for variables interfacing with the 'spectral_wave_data' API you need to set:
    kind_swd_c = c_double
else
    kind_swd_c = c_float
in 'kind_values.f90' before compiling the fortran library.

If the C-program apply 'float' for variables interfacing with the 'spectral_wave_data' API you also need to 
define the SWD_API_FLOAT macro before including 'SpectralWaveData.h'.

---------------------------------------------------------------------------------------------
Do the following steps on Windows:
---------------------------------------------------------------------------------------------
1) cd Cmake
2) Adjust CMakeBuild_Win64.bat to reffer to your actual version of Visual Studio
3) Run CMakeBuild_Win64.bat for creating the Visual Studio Solution
4) cd Build_Win64
5) Open swd_cpp_app_example.sln with Visual Studio
6) Build -> Build Solution
7) Check that the solution build without error messages
8) Set 'swd_cpp_app_example' project as startup project before you can run it in the GUI.
9) Run cleanup.bat when you have tested the program

---------------------------------------------------------------------------------------------
Do the following steps on Linux:
---------------------------------------------------------------------------------------------

1) cd Cmake
2) Run CMakeBuild_Linux.bash for creating the Makefile
3) cd Build_Linux
4) make


   