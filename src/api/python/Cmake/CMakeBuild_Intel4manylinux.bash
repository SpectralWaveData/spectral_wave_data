rm -rf ./Build_4manylinux
mkdir ./Build_4manylinux
cd ./Build_4manylinux
cmake -DCMAKE_Fortran_COMPILER=ifort -DSTATIC=OFF ..
make
