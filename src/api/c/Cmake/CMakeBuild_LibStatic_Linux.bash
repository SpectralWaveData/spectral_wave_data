rm -rf Build_Linux_LibStatic
mkdir Build_Linux_LibStatic
cd Build_Linux_LibStatic
cp ../CMakeLists_LibStatic.txt CMakeLists.txt
cmake -DCMAKE_Fortran_COMPILER=gfortran-10 -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_BUILD_TYPE=Release .
make
cd ..
python make_zip_Linux_LibStatic.py
