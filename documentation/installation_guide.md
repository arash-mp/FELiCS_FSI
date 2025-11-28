# Installation of FELiCS2.0

This installation guide provides step-by-step instructions for setting up the required software and dependencies for running FELiCS. It covers prerequisites, installation commands, and configuration steps to ensure a smooth setup process.

## Requirements

FELiCS2.0 is designed to run on UNIX-based systems, including native Linux (recommended distributions: Ubuntu, CentOS, Fedora) and macOS. If you have linux or Mac System, you can jump to section [Setup](#setup).

### For Windows Users
For Windows users, FELiCS2.0 can be used via the Windows Subsystem for Linux (WSL).

#### Step 1: Install WSL on your system: 
Follow the official [Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install) to set up WSL. **We strongly recommend installing WSL version 2 for improved performance and compatibility.** 

#### Step 2: Install a Linux distribution in your WSL
After setting up WSL select and install a Linux distribution. This is also coverend in the [Microsoft WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install). Theoretically, FELiCS should work with every recent Ubuntu version. However we tested it with the following: 
- Ubuntu 20.04
- Ubuntu 22.04
- Ubuntu 22.04.5 LTS


## Setup

### Step 1: Clone FELiCS Code
First, navigate to the location where you want to store the FELiCS folder and clone the git repository [FELiCS2.0](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0) by executing the following command line: 
```bash
git clone https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0
```
### Step 2: Install conda
The recommended way of installation uses [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for package management. Please make sure it is available on your machine or install it if not.


## Package Installation

### Step 1: Install conda Environment
In the installation directory, we provide different _.yml_ files that contain all required packages. The latest version that works on mutliple systems is 
- [felics_v2.4_env.yml](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/installation/felics_v2.4_env.yml)

The older versions of yml files can be found in the folder [yml_old](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/installation/yml_old) with the corresponding version in the filename. 


```{admonition} Optional
In case you want to give the environment a different name, you can change the first line of the _.yml_ file accordingly. name: `felics  ->  name: new_name`
```
As an example, if the downloaded _.yml_ file "_felics_v2.4_env.yml_" is in the "_yml_" folder of the git repository (located in your home folder), the command would be 
```bash
conda env create -f ~/felics2.0/installation/felics_v2.4_env.yml -y
conda env create -f ~/felics2.0/installation/felics_v2.4_env.yml -y
```

```{admonition} Optional
if you want to review the packages being installed, omit the `-y`
```

```{attention} 
if the installation fails due to "No space left on device", you can specify an alternative installation directory using the `--prefix` option:
`conda env create --prefix <path_to_more_space> -f ~/felics2.0/installation/felics2.0_env.yml -y`
```

### Step 2: Set FELiCS alias

An alias is a shortcut command that you define in your shell configuration (such as `.bashrc`). It allows you to run FELiCS by simply typing `FELiCS` in your terminal, instead of specifying the full path and command each time. This makes it easier and faster to launch FELiCS from any location in your terminal.

```{note}
For Mac users, the ```bashrc``` file is called ```zshrc```.
```

To create an alias, you have two options:

**Option 1: Automatically add the FELiCS alias to your `.bashrc`**

Now, go to the top level of the cloned git repository. The top level is the 'felics2.0' folder and it contains for example the 'src' folder and the 'REAMDE.md' file. 

Run the following command block in your terminal to add a FELiCS alias to your _~/.bashrc_:
```bash
export FELICS_DIRECTORY=`pwd`
if [ -f "$FELICS_DIRECTORY/src/main.py" ]; then
	echo -e "FELiCS_PATH=\"${FELICS_DIRECTORY}\" 
	FELiCS() {
	    export OMP_NUM_THREADS=2;
	    python \$FELiCS_PATH/src/main.py \"\$@\";
	    python \$FELiCS_PATH/src/main.py \"\$@\";
	    unset OMP_NUM_THREADS
	}" >> ~/.bashrc
else
	echo -e "\033[0;31m Error: this does not seem to be a FELiCS directory.\033[0m"
fi
```

```{caution}
if you receive the error "`this does not seem to be a FELiCS directory`", please make sure that you are in the felics folder. The top level is the `felics2.0` folder and it contains for example the `src` folder and the `REAMDE.md` file. 
```

After adding content to your `~/.bashrc`, you need to run it. To do so, run ```source ~/.bashrc```

**Option 2: Manually add the lines to your .bashrc**

Instead of running the above bash command block, you can also add the FELiCS alias manually by inserting the following paragraph into your _~/.bashrc_

```bash
FELiCS_PATH="<path_to_felics_repository>" 
FELiCS() {
    export OMP_NUM_THREADS=2;
    python $FELiCS_PATH/src/main.py "$@" ;
    unset OMP_NUM_THREADS
}
```
After adding content to your `~/.bashrc`, you need to run it. To do so, run ```source ~/.bashrc```


### Step 3: Verify Installation

After completing the installation steps, open a new terminal and activate your conda environment:
```bash
conda activate felics
```
If you used a different environment name, activate it accordingly:
```bash
conda activate <your_environment_name>
```

Once the environment is active, test the FELiCS alias:
```bash
FELiCS
```
If successful, the FELiCS logo is displayed:
```bash
(         (               (
)\ )      ) )       (    )\ )
(()/(  (  (()/( (    )\   (()/(
/(_)) )\  /(_)))\  (((_)  /(_))
(_)_)((_) (_)) ((_) )\___ (_))
| __|| __|| |   (_)((/ __|/ __|
| _| | _| | |__ | | | (__ \__ \
|_|  |___||____||_|  \___||___/
```


For a more thorough check, run a tutorial case:
```bash
FELiCS -f $FELiCS_PATH/tutorials/modal_analysis_tutorial/modal.json
```
This example should complete in under a minute. If successful, you will see the message: `Finished FELiCS run.`


## Installation of FELiCS as a package (optional)
To use FELiCS as a regular python package for scripting, it can be installed in the conda environment using pip and then easily be imported. 

```{admonition} Optional: Install FELiCS as a Python Package
If you do not plan to use FELiCS as a Python package within your own scripts, you can skip this step.
```
### Step 1: Install FELiCS as a python package

To install FELiCS as a python package, run the following commands:

```bash
cd </path/to/felics>/installation
conda activate <felics-environemnt>
pip install -e .
```

```{note}
The flag `-e` indicates an installation in editable mode. This allows the FELiCS repository to be connected to the conda environment and all changes are reflected instantaneously in the environment.
```

### Step 2: Verify installation

To test the installation, start python in the FELiCS conda environment

```bash
conda activate <felics-environemnt>
python import FELiCS
```
and try to import FELiCS. If successful, the FELiCS logo is displayed:
```python import FELiCS
(         (               (
)\ )      ) )       (    )\ )
(()/(  (  (()/( (    )\   (()/(
/(_)) )\  /(_)))\  (((_)  /(_))
(_)_)((_) (_)) ((_) )\___ (_))
| __|| __|| |   (_)((/ __|/ __|
| _| | _| | |__ | | | (__ \__ \
|_|  |___||____||_|  \___||___/
```

```{Caution}
The package name is case sensitive (`FELiCS`). In the conda environment, it may appear in lowercase, but you must use the correct capitalization when importing.
```

If all these steps complete without errors, your FELiCS installation is ready to use.


If installation using pip is not preferred, FELiCS can also be added to the `$PYTHONPATH` in the `.bashrc`. To do this add the line
```bash
export PYTHONPATH=$PYTHONPATH:$FELiCS_PATH
```
to the `.bashrc`.

## Uninstalling FELiCS
Because FELiCS is installed in editable mode the standard unistallation procedure does not work. To uninstall the folder `$FELiCS_PATH/src/FELiCS.egg-info` must be removed. In the folder of the conda environment (`$CONDA_PATH/envs/<felics-environment>`) search for felics and delete all files containing that name. Like this the package is fully removed from the environemnt.

After changing your `~/.bashrc` you need to restart the terminal or ssh connection to make the changes take effect. 

## Additional Software

```{admonition} Optional
The following software packages are not required to run the FELiCS code, but it improves the tutorial experience and will be usefull for most users.
```
For enhanced pre- and postprocessing, and to complete all tutorial tasks, the following software is recommended:

- **[GMSH](https://gmsh.info/) with graphical user interface:**  
	A mesh generation tool with a graphical user interface.  
	- *Windows users:* You can install GMSH directly on Windows, even if FELiCS is running in WSL.

- **[ParaView](https://www.paraview.org/):**  
	A powerful visualization tool for viewing data on structured and unstructered grids.  
	- *Recommended for tutorials:* Some tutorial steps require ParaView for full functionality.
	- *Alternatives:* You may also post-process data using `matplotlib` or your own Python scripts.

