rmdir /s /q Build_Win64
mkdir Build_Win64
cd Build_Win64
rem cmake -G "Visual Studio 15 2017 Win64" ..
cmake -G "Visual Studio 16 2019" -A x64 ..
pause
