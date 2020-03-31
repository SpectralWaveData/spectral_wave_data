rm -rf Build_Linux
mkdir Build_Linux
cd Build_Linux
cmake -DCMAKE_Fortran_COMPILER=gfortran-9 -DCMAKE_C_COMPILER=gcc-9 -DSTATIC=ON -DCMAKE_BUILD_TYPE=Release ..
make
