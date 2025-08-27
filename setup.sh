#!bin/bash
# for use on a debian-12 root machine
apt update
apt install cuda-toolkit

wget https://developer.download.nvidia.com/compute/cuda/13.0.0/local_installers/cuda_13.0.0_580.65.06_linux.run
sh cuda_13.0.0_580.65.06_linux.run


wget ftp://ftp.gromacs.org/gromacs/gromacs-2025.2.tar.gz

tar xfz gromacs-2025.2.tar.gz 
cd gromacs-2025.2

mkdir build 
cd build

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda

make
make check
make install
source /usr/local/gromacs/bin/GMXRC