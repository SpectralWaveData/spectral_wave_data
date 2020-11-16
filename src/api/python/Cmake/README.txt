For building on Windows with static linking of Intel libraries into SpectralWaveData.dll:
=========================================================================================
This option removes the need for Intel Redistributable libraries.

1) Run CmakeBuild_Win64.bat
2) Open "Build_Win64\SpectralWaveData.sln in Visual Studio
3) Change from debug to release configuration
4) Set SpectralWaveData as current project
5) Open Projects -> Properties -> Fortran -> Libraries -> Runtime Library
   Select value "Multithreaded"   Reference for more info can be found in this link:
   https://community.intel.com/t5/Intel-Fortran-Compiler/DLL-generation-and-static-linking-to-intel-runtime-libraries/td-p/1183187
6) Build project

For building on Windows with dynamic linking of Intel libraries into SpectralWaveData.dll:
=========================================================================================
This option may require installation of Intel Redistributable libraries or Visual studio.

Same steps as above, except skip step 5.

For building on Linux with static linking of Intel libraries into SpectralWaveData.so
=====================================================================================
The Cmake files just works :-)
