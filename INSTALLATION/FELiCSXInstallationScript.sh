sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update

# dolfinx with real numbers:
sudo apt install fenicsx

# add complex number support
sudo apt install python3-dolfinx-complex

# add to bashrc:
export PETSC_DIR=/usr/lib/petscdir/petsc-complex
export SLEPC_DIR=/usr/lib/slepcdir/slepc-complex

# add h5py
sudo apt-get install python3-h5py python3-h5py-mpi
