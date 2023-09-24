if "%1"== ""   set generator="Visual Studio 17 2022" -A x64
if "%1"== "17" set generator="Visual Studio 17 2022" -A x64
if "%1"== "16" set generator="Visual Studio 16 2019" -A x64
if "%1"== "15" set generator="Visual Studio 15 2017 Win64"

rmdir /s /q Build_Win64_LibStatic
mkdir Build_Win64_LibStatic
cd Build_Win64_LibStatic
copy ..\CMakeLists_LibStatic.txt CMakeLists.txt
cmake -G %generator% .
pause
