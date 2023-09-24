rm -rf Build_Linux
mkdir Build_Linux
cd Build_Linux
cmake -DCMAKE_Fortran_COMPILER=gfortran-11 -DCMAKE_BUILD_TYPE=Release ..
make
