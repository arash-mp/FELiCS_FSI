# Installation of FELiCS2.0

## Requirements

FELiCS2.0 is designed to run on UNIX-based systems, including native Linux (recommended distributions: Ubuntu, CentOS, Fedora) and macOS.

For Windows users, FELiCS2.0 can be used via the Windows Subsystem for Linux (WSL). Follow the official [Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install) to set up WSL. **We strongly recommend installing WSL version 2 for improved performance and compatibility.** After setting up WSL select and install a Linux distribution. We recommend Ubuntu (latest LTS version).

## Setup
First, navigate to the location where you want to store the FELiCS folder and clone the git repository [FELiCS2.0](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0) by executing the following command line: 
```bash
git clone https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0
```

The recommended way of installation uses [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for package management. Please make sure it is available on your machine or install it if not.

After cloning the repository and installing conda, make sure you take the following three steps:
1. [Install Necessary Packages](#package-installation)
2. [Add FELiCS alias](#felics-alias)
3. [Verify Installtion](#Verify-Installation)


## Package Installation
### Conda Environment
In the INSTALLATION directory, we provide different _.yml_ files that contain all required packages. The latest version that works on mutliple systems is 
- [felics_v2.4_env.yml](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/INSTALLATION/felics_v2.4_env.yml)

The older versions of yml files can be found in the folder [yml_old](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/INSTALLATION/yml_old) with the corresponding version in the filename. 

>**Optional:** in case you want to give the environment a different name, you can change the first line of the _.yml_ file accordingly
>```
>name: felics2025_dolfin9  ->  name: new_name
>```

To create the conda environment, enter the following command in the terminal:
```bash
conda env create -f <path_to_your_yml_file> -y
```
As an example, if the downloaded _.yml_ file "_felics_v2.4_env.yml_" is in the "_yml_" folder of the git repository (located in your home folder), the command would be 
```bash
conda env create -f ~/felics2.0/INSTALLATION/felics_v2.4_env.yml -y
```
>**Optional:** if you want to review the packages being installed, omit the `-y`

>**Remark:** if the installation fails due to "No space left on device", you can specify an alternative installation directory using the `--prefix` option:
>```bash
>conda env create --prefix <path_to_more_space> -f ~/felics2.0/INSTALLATION/felics2.0_env.yml -y
>```

## FELiCS alias
Now, go to the top level of the cloned git repository and run this command block in your terminal to add a FELiCS alias to your _~/.bashrc_:
```bash
export FELICS_DIRECTORY=`pwd`

if [ -f "$FELICS_DIRECTORY/src/main.py" ]; then
	echo -e "FELiCS_PATH=\"${FELICS_DIRECTORY}\" 
	FELiCS() {
	    export OMP_NUM_THREADS=2;
	    python \$FELiCS_PATH/src/main.py \"\$@\";
	    unset OMP_NUM_THREADS
	}" >> ~/.bashrc
else
	echo -e "\033[0;31m Error: this does not seem to be a FELiCS directory.\033[0m"
fi
```
>**Remark:** if you receive the error "`this does not seem to be a FELiCS directory`", please make sure that you are in the folder containing the FELiCS _main.py_.

>**Optional:** instead of running the above bash command block, you can also add the FELiCS alias manually by inserting the following paragraph into your _~/.bashrc_
```bash
FELiCS_PATH="<path_to_felics_repository>" 
FELiCS() {
    export OMP_NUM_THREADS=2;
    python $FELiCS_PATH/main.py "$@" ;
    unset OMP_NUM_THREADS
}
```
Do not forget to source it with ```source ~/.bashrc```

## Installation of FELiCS as a package
To use FELiCS as a regular python package, it can be installed in the conda environment using pip and then easily be imported. To install FELiCS as a python package, run the following commands:

```bash
cd /path/to/felics
conda activate <felics-environemnt>
pip install -e .
```

```{note}
The flag `-e` indicates an installation in editable mode. This allows the FELiCS repository to be connected to the conda environment and all changes are reflected instantaneously in the environment.
```

To test the installation, start python in the FELiCS conda environment

```bash
conda activate <felics-environemnt>
python
```
and try to import FELiCS. If successful, the FELiCS logo is displayed:
```python
import FELiCS
(         (               (
)\ )      ) )       (    )\ )
(()/(  (  (()/( (    )\   (()/(
/(_)) )\  /(_)))\  (((_)  /(_))
(_)_)((_) (_)) ((_) )\___ (_))
| __|| __|| |   (_)((/ __|/ __|
| _| | _| | |__ | | | (__ \__ \
|_|  |___||____||_|  \___||___/
```

If installation using pip is not preferred, FELiCS can also be added to the `$PYTHONPATH` in the `.bashrc`. To do this add the line
```bash
export PYTHONPATH=$PYTHONPATH:$FELiCS_PATH
```
to the `.bashrc`.

### Uninstall
Because FELiCS is installed in editable mode the standard unistallation procedure does not work. To uninstall the folder `$FELiCS_PATH/src/FELiCS.egg-info` must be removed. In the folder of the conda environment (`$CONDA_PATH/envs/<felics-environment>`) search for felics and delete all files containing that name. Like this the package is fully removed from the environemnt.

After changing your `~/.bashrc` you need to restart the terminal or ssh connection to make the changes take effect. 

## Verify Installation

After completing the installation steps, open a new terminal and activate your conda environment:
```bash
conda activate felics2025_dolfin9
```
If you used a different environment name, activate it accordingly:
```bash
conda activate <your_environment_name>
```

Once the environment is active, test the FELiCS alias:
```bash
FELiCS
```
You should see the FELiCS logo displayed.

For a more thorough check, run a tutorial case:
```bash
FELiCS -f $FELiCS_PATH/TUTORIALS/modal_analysis_tutorial/modal.json
```
This example should complete in under a minute. If successful, you will see the message: `Finished FELiCS run.`

To verify the FELiCS package installation, run:
```bash
python -c "import FELiCS"
```
> **Note:** The package name is case sensitive (`FELiCS`). In the conda environment, it may appear in lowercase, but you must use the correct capitalization when importing.

If all these steps complete without errors, your FELiCS installation is ready to use.


