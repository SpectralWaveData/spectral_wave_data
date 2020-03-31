The kind_values.f90 file is important to ensure that the numerical precision
(float or double) applied in the C-program match the fortran code.

If the C-program apply 'double' for variables interfacing with the 'spectral_wave_data' API you need to set:
    kind_swd_c = c_double
else
    kind_swd_c = c_float
in 'kind_values.f90' before compiling the fortran library.

If the C-program apply 'float' for variables interfacing with the 'spectral_wave_data' API you also need to 
define the SWD_API_FLOAT macro before including 'spectral_wave_data.h'.

---------------------------------------------------------------------------------------------
Do the following steps on Windows:
---------------------------------------------------------------------------------------------
1) cd Cmake
2) Run CMakeBuild_Win64.bat for creating the Visual Studio Solution
3) cd Build_Win64
4) Open swd_c_app_example.sln with Visual Studio
5) Build -> Build Solution
6) Check that the solution build without error messages
7) Set 'swd_c_app_example' project as startup project before you can run it in the GUI.
8) Run cleanup.bat when you have tested the program

---------------------------------------------------------------------------------------------
Do the following steps on Linux:
---------------------------------------------------------------------------------------------

1) cd Cmake
2) Run CMakeBuild_Linux.bash for creating the Makefile
3) cd Build_Linux
4) make



