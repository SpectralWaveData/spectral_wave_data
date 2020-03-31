rmdir /s /q Build_Win64_LibStatic
mkdir Build_Win64_LibStatic
cd Build_Win64_LibStatic
copy ..\CMakeLists_LibStatic.txt CMakeLists.txt
cmake -G "Visual Studio 15 2017 Win64" .
pause
