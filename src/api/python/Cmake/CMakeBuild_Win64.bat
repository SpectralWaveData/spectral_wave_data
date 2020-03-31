rmdir /s /q Build_Win64
mkdir Build_Win64
cd Build_Win64
cmake -G "Visual Studio 15 2017 Win64" -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE ..
pause
