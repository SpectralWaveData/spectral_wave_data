rm -rf  build
rm -rf  dist
rm -rf  spectral_wave_data
rm -rf  spectral_wave_data.egg-info
rm -rf  MANIFEST.in

mkdir spectral_wave_data
cp -r ../spectral_wave_data/* spectral_wave_data/.
cp -r ../Cmake/Build_Linux/* spectral_wave_data/.
cp -r ../../c/spectral_wave_data.h spectral_wave_data/.
cp -r ../../../../icons/SWD_logo.ico spectral_wave_data/.

cp LICENSE.txt spectral_wave_data/.
cp README.txt spectral_wave_data/.

python setup.py bdist_wheel
cp -r dist/* ../whl/.

echo ""
echo "You may want to run python -m auditwheel repair whl/*"
echo "to create a manylinux wheel"
