# Input files

The main input files for FELiCS are:
- a configuration file (`.json`)
- a mesh file (`.msh`)
- a boundary conditions file (`.json`)
- a mean flow file (`.fel` or `MeanFlow.h5`)

Additional input files may be required depending on the simulation setup.

The different input files are described in the following sections:

| File | Description |
| - | - |
| [Configuration file](param_files.md) | Main simulation settings and solver parameters |
| [Mesh file](gmsh_files.md) | Computational mesh in Gmsh format |
| [Mean flow file](mean_files.md) | Mean flow solution used by FELiCS |
| [Boundary file](bc_files.md) | Boundary condition definitions |
| [Mixture file](mix_files.md) | Optional mixture properties for reacting flows |

```{toctree}
:maxdepth: 2

param_files
gmsh_files
mean_files
bc_files
mix_files
```