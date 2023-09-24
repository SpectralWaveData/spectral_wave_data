rmdir /s /q build
rmdir /s /q dist
rmdir /s /q spectral_wave_data
rmdir /s /q spectral_wave_data.egg-info
del /q MANIFEST.in

mkdir spectral_wave_data
xcopy /e ..\spectral_wave_data\* spectral_wave_data\*
copy ..\Cmake\Build_Win64\Release\* spectral_wave_data\*
copy ..\..\c\spectral_wave_data.h spectral_wave_data\*
copy ..\..\..\..\icons\SWD_logo.ico spectral_wave_data\*

copy LICENSE.txt spectral_wave_data\*
copy README.txt spectral_wave_data\*

python setup.py bdist_wheel
copy dist\* ..\whl\*

pause
