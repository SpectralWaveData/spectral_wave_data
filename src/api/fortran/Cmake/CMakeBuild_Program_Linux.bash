rm -r -f Build_Linux
mkdir Build_Linux
cd Build_Linux
cmake -DCMAKE_Fortran_COMPILER=gfortran-9 -LIBTYPE=STATIC ..
make