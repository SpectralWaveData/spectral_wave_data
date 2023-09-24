import os
import shutil

"""
Copy all relevant files for building manylinux distro to a subfolder "files4manylinux".
"""

zip_name = "files4manylinux"

zip_dir = zip_name
package_dir = os.path.join(zip_dir, "spectral_wave_data")

if os.path.isdir(zip_dir):
    shutil.rmtree(zip_dir)

os.mkdir(zip_dir)

shutil.copytree(os.path.join("..", "spectral_wave_data"), package_dir)
shutil.copy(os.path.join("..", "Cmake", "Build_4manylinux", "libSpectralWaveData.so"), package_dir)

shutil.copy("setup.py", zip_dir)
shutil.copy("readme_pypi.md", zip_dir)
shutil.copy("VERSION.dat", zip_dir)
shutil.copy("LICENSE.txt", package_dir)
shutil.copy("README.txt", package_dir)
shutil.copy(os.path.join("..", "..", "c", "spectral_wave_data.h"), package_dir)
shutil.copy(os.path.join("..", "..", "..", "..", "icons", "SWD_logo.ico"), package_dir)

cp ../../c/spectral_wave_data.h spectral_wave_data/.
cp ../../../../icons/SWD_logo.ico spectral_wave_data/.


shutil.make_archive(zip_name, 'gztar', zip_dir)
#shutil.make_archive(zip_name, 'zip', zip_dir)

"""
On your manylinux distro:
Unpack files4manylinux.tar.gz in a temporary build directory
cd to_build_directory_containing_setup.py
python setup.py bdist_wheel
python -m auditwheel repair dist/*.whl

The previous step generates the final whl to be uploaded on test.pypi.org. 
On any linux system try:
    "pip install --index-url https://test.pypi.org/simple/ spectral_wave_data"
If the package work you can upload the whl file to the real pypi.org
site for official distribution.
"""
