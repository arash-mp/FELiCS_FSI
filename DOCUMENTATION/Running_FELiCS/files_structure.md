# Files needed to run FELiCS

## Mesh
Your case folder must containa mesh in the **.msh** format. 
This can be generated for instance with GMSH. Make sure the .msh file is saved in *Version 2 ASCII* format.
 An example for the mesh generation is provided in [Tutorial 1: Helical Sphere Wake Instability](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/main/DOCUMENTATION/_build/Tutorials/tutorial_1).

## Mean flow

##

Once the installation is completed, all the input parameters will only be loaded from the three following **.json** files:
- `settings.json`: contains most of the information about the FELiCS run.
- `boundaries.json`: con*string*tains specific information about the boundary conditions.
- `mixture.json`: contains specific information about the physical properties of the fluid considered.
