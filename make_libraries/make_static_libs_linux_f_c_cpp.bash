dir_libs=$(pwd)/libs_linux_f_c_cpp

dir_f=$(pwd)/../src/api/fortran/Cmake
dir_c=$(pwd)/../src/api/c/Cmake
dir_cpp=$(pwd)/../src/api/cpp/Cmake

cd $dir_f
bash CMakeBuild_LibStatic_Linux.bash
cp swd_lib_static_f_linux_*.tar.gz $dir_libs/.

cd $dir_c
bash CMakeBuild_LibStatic_Linux.bash
cp swd_lib_static_c_linux_*.tar.gz $dir_libs/.

cd $dir_cpp
bash CMakeBuild_LibStatic_Linux.bash
cp swd_lib_static_cpp_linux_*.tar.gz $dir_libs/.
