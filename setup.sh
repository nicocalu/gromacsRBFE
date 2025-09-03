#!/bin/bash
# to run on a debian-12 root with at least 1 CUDA capable gpu
set -euo pipefail

if [[ $EUID -ne 0 ]]; then
  echo "Run as root."
  exit 1
fi

apt-get update
apt-get install -y build-essential dkms linux-headers-$(uname -r) git wget curl pv pigz python3.11-venv cmake

cd /tmp

# Install NVIDIA driver from CUDA repo
wget -q https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/cuda-keyring_1.1-1_all.deb
dpkg -i cuda-keyring_1.1-1_all.deb
apt-get update
apt-get install -y cuda-drivers

echo "NVIDIA driver installed. Reboot required to load the kernel module."
echo "After reboot, run: /tmp/continue_gromacs_setup.sh"
cat >/tmp/continue_gromacs_setup.sh <<'EOS'
#!/bin/bash
set -euo pipefail

# Verify driver
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi not found; ensure driver installed and Secure Boot not blocking."
  exit 1
fi
nvidia-smi || { echo "Driver not loaded. Check Secure Boot or logs."; exit 1; }

cd /tmp

# NVIDIA HPC SDK
wget https://developer.download.nvidia.com/hpc-sdk/25.7/nvhpc_2025_257_Linux_x86_64_cuda_12.9.tar.gz
pv nvhpc_2025_257_Linux_x86_64_cuda_12.9.tar.gz | pigz -dc -p8 | tar -xf -

./nvhpc_2025_257_Linux_x86_64_cuda_12.9/install

# Environment 
NVHPC=/tmp/nvidia/hpc_sdk

cat >/etc/profile.d/hpcsdk.sh <<'EOF'
# NVIDIA HPC SDK environment
export NVHPC=/tmp/nvidia/hpc_sdk
export HPCSDK="$NVHPC"
export PATH="$NVHPC/Linux_x86_64/25.7/compilers/bin:$NVHPC/Linux_x86_64/25.7/cuda/12.9/bin:$PATH"
# Ensure defined under 'set -u'
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$NVHPC/Linux_x86_64/25.7/compilers/lib:$NVHPC/Linux_x86_64/25.7/cuda/12.9/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CUDA_HOME="$NVHPC/Linux_x86_64/25.7/cuda/12.9"
EOF
chmod +x /etc/profile.d/hpcsdk.sh
# Load env in current shell
source /etc/profile.d/hpcsdk.sh

# GROMACS
wget -q ftp://ftp.gromacs.org/gromacs/gromacs-2025.2.tar.gz
tar -xzf gromacs-2025.2.tar.gz
cd gromacs-2025.2
mkdir -p build && cd build

cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=CUDA \
  -DGMX_SIMD=AVX2_256 \
  -DCMAKE_INSTALL_PREFIX=/tmp/gromacs \
  -DGMX_NVSHMEM=ON \
  -DNVSHMEM_ROOT=$NVHPC/Linux_x86_64/25.7/comm_libs/nvshmem

make -j"$(nproc)"
make check -j"$(nproc)"
make install

source /tmp/gromacs/bin/GMXRC
echo 'source /tmp/gromacs/bin/GMXRC' >> /etc/profile.d/gromacs.sh
chmod +x /etc/profile.d/gromacs.sh
source /tmp/gromacs/bin/GMXRC

echo "GROMACS installed. Open a new shell or 'source /etc/profile' to load env."

cd /tmp
git clone https://github.com/nicocalu/gromacsRBFE
cd gromacsRBFE

python3 -m venv grbfe
source grbfe/bin/activate
python3 -m pip install -r requirements.txt

EOS
chmod +x /tmp/continue_gromacs_setup.sh

echo "Reboot now, then run: sudo /tmp/continue_gromacs_setup.sh"
# gmx pdb2gmx -f compound.pdb
# gmx editconf -f conf.gro -bt dodecahedron -d 1.0 -o boxed.gro
# gmx solvate -cp boxed.gro -cs spc216 -p topol.top -o solvated.pdb
# gmx genion -s -nn 1 -o ions.pdb -nname CL -pname NA -p