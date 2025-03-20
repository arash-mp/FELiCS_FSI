#!/bin/bash
# This file installs FELiCS in on the cluster at HFI. Please make sure you have homebrew installed
# Furthermore, you need access to the FELiCS gitlab repository. To get it please
# email t.kaiser@tu-berin.de. Once you have access to the online repository add the
# ssh-key of your mashine to your gitlab key-list (see
# https://docs.gitlab.com/ee/ssh/README.html#adding-an-ssh-key-to-your-gitlab-account ).
# Then source this file from the directory you want to install FELiCS into. You may 
# delete the install file afterwards. 
module load python3
echo "Create FELiCS environement in anaconda"
conda create -n FELiCS
source activate FELiCS
conda install scipy
conda install matplotlib
conda install spyder
conda install colorama
conda install h5py
conda config --set channel_priority strict
conda install -c conda-forge meshio=3.3.0
conda install -c conda-forge fenics
pip install pyvtk

git clone git@gitlab.tubit.tu-berlin.de:thomaskaiser/FELiCS.git

echo "Defining FELiCS command by adding it to \~.bashrc"
echo "export PYTHONPATH=\$PYTHONPATH:’"$PWD"/FELiCS’" >> ~/.bashrc
echo 'FELiCS() {' >> ~/.bashrc
echo 'export OMP_NUM_THREADS=2;' >> ~/.bashrc
echo 'python '$PWD'/FELiCS/main.py ”$@”;' >> ~/.bashrc
echo 'unset OMP_NUM_THREADS' >> ~/.bashrc
echo '}' >> ~/.bashrc
source ~/.bashrc
