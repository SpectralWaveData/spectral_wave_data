rmdir /s /q Build_Win64_LibStatic
mkdir Build_Win64_LibStatic
cd Build_Win64_LibStatic
copy ..\CMakeLists_LibStatic.txt CMakeLists.txt
rem cmake -G "Visual Studio 15 2017 Win64" .
cmake -G "Visual Studio 16 2019" -A x64 .
pause
