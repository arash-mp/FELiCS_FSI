# Files needed to run FELiCS

## Case folder
In the beginning of every set up of a FELiCS case, a case folder should be created. In this folder the setting file(s) as well as the output directory or directories for the results will be located. Depending on the case it might also include the input files.
Here is an example of a FELiCS case directory content:

<pre>
myFELiCScase/
├── settings.json
├── boundary.json
└── Base_flow/
    ├── base_flow.fel
└── Mesh/
    ├── mesh.msh
└── Outputs/ </pre>

## Mesh
Your case folder must containa mesh in the **.msh** format. 
This can be generated for instance with GMSH. Make sure the .msh file is saved in *Version 2 ASCII* format. FeliCS only processes triangular cells, make sur to not "Recombine" you mesh with GMSH. Physical entities must be defined in the `.msh `file in order to impose boundary conditions.
 An example for the mesh generation is provided in [Tutorial 1: Helical Sphere Wake Instability](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/documentation/DOCUMENTATION/_build/Tutorials/tutorial_1).

## Base flow
The base flow can be obtained by various sources depending on the case: Numerical simulations, such as RANS, URANS, LES, DNS, as well as experimental results or analytical models. A RANS mean field can also be computed using the finite element Newton solver `FlowSolver.py` integrateed in FELiCS (see [Tutorial 1: Helical Sphere Wake Instability](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/documentation/DOCUMENTATION/_build/Tutorials/tutorial_1)).

The mean flow information (velocities, pressure, viscosity, forcing and response domains, ...) are encapsulated in a **.fel** file. 
The construction of this file is presented in [fel file in FELiCS](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/documentation/DOCUMENTATION/_build/Running_FELiCS/fel_file).

## Analysis parameters

All the input parameters are loaded from the following **.json** files:
- `settings.json`: contains the analysis settings and files directory.
- `boundaries.json`: contains the boundary conditions.

## Output directory
Prepare an ouptuts folder where the FELiCS result files will be outputed.
