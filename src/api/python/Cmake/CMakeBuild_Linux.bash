rm -rf Build_Linux
mkdir Build_Linux
cd Build_Linux
cmake -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DSTATIC=OFF ..
make
