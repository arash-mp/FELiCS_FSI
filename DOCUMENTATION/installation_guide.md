# Installation of FELiCS2.0
## Setup
First, navigate to the location where you want to store the FELiCS folder and clone the git repository [FELiCS2.0](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0) by executing the following command line: 
```bash
git clone https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0
```

The recommended way of installation uses [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for package management. Please make sure it is available on your machine or install it if not.

After cloning the repository and installing conda, make sure you take the following two steps:
1. [install necessary packages](#package-installation)
2. [add FELiCS alias](#felics-alias)


## Package Installation
### Conda Environment
In the yml directory, we provide different _.yml_ files that contain all required packages. The latest version that worked on mutliple systems is 
- [felics2.0_env.yml](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/main/yml/felics2.0_env.yml)

Older versions of yml files can be found in the folder [yml/old_versions](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/tree/main/yml/old_versions?ref_type=heads)

>**Optional:** in case you want to give the environment a different name, you can change the first line of the _.yml_ file accordingly
>```
>name: felics  ->  name: new_name
>```

To create the conda environment, enter the following command in the terminal:
```bash
conda env create -f <path_to_your_yml_file> -y
```
As an example, if the downloaded _.yml_ file "_felics2.0_env.yml_" is in the "_yml_" folder of the git repository (located in your home folder), the command would be 
```bash
conda env create -f ~/felics2.0/yml/felics2.0_env.yml -y
```
>**Optional:** if you want to review the packages being installed, omit the `-y`

>**Remark:** if the installation fails due to "No space left on device", you can specify an alternative installation directory using the `--prefix` option:
>```bash
>conda env create --prefix <path_to_more_space> -f ~/felics2.0/yml/felics2.0_env.yml -y
>```

## FELiCS alias
Now, go to the top level of the cloned git repository and run this command block in your terminal to add a FELiCS alias to your _~/.bashrc_:
```bash
export FELICS_DIRECTORY=`pwd`

if [ -f "$FELICS_DIRECTORY/main.py" ]; then
	echo -e "FELiCS_PATH=\"${FELICS_DIRECTORY}\" 
	FELiCS() {
	    export OMP_NUM_THREADS=2;
	    python \$FELiCS_PATH/main.py \"\$@\";
	    unset OMP_NUM_THREADS
	}" >> ~/.bashrc
else
	echo -e "\033[0;31m Error: this does not seem to be a FELiCS directory.\033[0m"
fi
```
>**Remark:** if you receive the error "`this does not seem to be a FELiCS directory`", please make sure that you are in the folder containing the FELiCS _main.py_.

>**Optional:** instead of running the above bash command block, you can also add the FELiCS alias manually by inserting the following paragraph into your _~/.bashrc_
>```
>FELiCS_PATH="<path_to_felics_repository>" 
>FELiCS() {
>    export OMP_NUM_THREADS=2;
>    python $FELiCS_PATH/main.py "$@" ;
>    unset OMP_NUM_THREADS
>}
>```

## Installation of FELiCS as a package
To use FELiCS as a package it can be installed in the conda environment using pip and then easily be imported. The flag `-e` indicates installation in editable mode. Like this the folder of the repo is connected to the environment and all changes are reflected instantaneously in the environment. Then run the following commands:
```
cd /path/to/felics
conda activate <felics-environemnt>
pip install -e .
```

If installation using pip is not preferred, FELiCS can also be added to the `$PYTHONPATH` in the `.bashrc`. To do this add the line
```
export PYTHONPATH=$PYTHONPATH:$FELiCS_PATH
```
to the `.bashrc`.

### Uninstall
Because FELiCS is installed in editable mode the standard unistallation procedure does not work. To uninstall the folder `$FELiCS_PATH/src/FELiCS.egg-info` must be removed. In the folder of the conda environment (`$CONDA_PATH/envs/<felics-environment>`` search for felics and delete all files containing that name. Like this the package is fully removed from the environemnt.

After changing your _~/.bashrc_ you need to restart the terminal or ssh connection to make the changes take effect. 


## Verify Installation
After finishing the installation steps you should be able to activate the new environment in your freshly opened terminal using
```bash
conda activate felics
```
>**Remark:** if you chose another name for your environment, use the same name for activating the environment:
>```bash
>conda activate <other_name>
>```

If that succeeds, try running FELiCS with 
```bash
FELiCS
```
>**Remark:** when working on c14 via ssh make sure to enable visual output with the `-X` option of ssh.

To verify the installation of FELiCS as a package run 
```bash
python -c "import FELiCS"
```
>**Remark:** Please note, that the package name is case sensitive, even though in the conda environment it is listed in lowercase letters. 

In case of an error, please refer to the chapter [common problems](#CommonIssues) for possible solutions.

The environment can be deactivated using 
```bash
conda deactivate felics
```

### Manual Package Installation
---------- **OUT OF DATE! SHOULD BE UPDATED TO NEWER PACKAGES** ------------
If the installation via _.yml_ files does not work, you can try installing the packages manually or [older yml files](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/tree/main/yml/old_versions?ref_type=heads). Since package compatibility can be tricky, it is recommended to stick to the following packages and versions which seem to be stable at the moment [August 2023].
```bash
conda create --name felics2.0 python=3.8
conda activate felics2.0
python3 -m pip cache purge
conda install numpy=1.23.3
pip install scipy==1.9.1 gmsh==4.10.5
conda install matplotlib=3.3.2=0 colorama=0.4.5 h5py=3.7.0
conda install -c conda-forge fenics-basix=0.5.0 fenics-dijitso=2019.1.0 fenics-dolfinx=0.5.1 fenics-fiat=2019.1.0
conda install -c conda-forge pyvista
```

