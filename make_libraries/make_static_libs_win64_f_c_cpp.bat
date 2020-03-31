rem You need to compile each library using Visual Studio before you can run this script to collect library zip files!

SET pwd=%CD%

SET dir_libs=%pwd%\libs_win64_f_c_cpp

SET dir_f=%pwd%\..\src\api\fortran\Cmake
SET dir_c=%pwd%\..\src\api\c\Cmake
SET dir_cpp=%pwd%\..\src\api\cpp\Cmake

cd %dir_f%
python make_zip_win64_LibStatic.py
copy /y swd_lib_static_f_win64_*.zip %dir_libs%\.

cd %dir_c%
python make_zip_win64_LibStatic.py
copy /y swd_lib_static_c_win64_*.zip %dir_libs%\.

cd %dir_cpp%
python make_zip_win64_LibStatic.py
copy /y swd_lib_static_cpp_win64_*.zip %dir_libs%\.

cd %pwd%
