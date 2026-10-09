# Mean-flow file

FELiCS imports mean-flow information using HDF5-based files.

Two types of mean-flow files are currently supported:
- `.fel` files containing mean-flow data defined on external CFD grids,
- `MeanFlow.h5` files previously exported by FELiCS on the computational mesh.

The mean-flow file used by a simulation is specified through [`FlowInput.MeanFlowFilePath`](parameters/FlowInput/MeanFlowFilePath.md).

## Workflows to import mean flows

FELiCS supports two approaches for importing mean-flow fields, depending on the value of the [`Case.needInterpolation`](parameters/Case/needInterpolation.md) parameter.

### 1. With interpolation

When `needInterpolation: true`, FELiCS assumes that the mean flow is defined on a grid different from the computational mesh.

The imported fields are then linearly interpolated onto the FEM mesh defined by [`MeshFilePath`](parameters/Case/MeshFilePath.md).

This is typically used with `.fel` files generated from:
- CFD simulations,
- experimental reconstructions,
- analytical models.

Example configuration:

```json
{
    "Case": {
        "needInterpolation": true,
        "MeshFilePath": "mesh.msh"
    },

    "FlowInput": {
        "MeanFlowFilePath": "meanflow.fel"
    }
}
```

For large meshes, particularly 3D cases, the interpolation step can become time-consuming.

After interpolation, FELiCS automatically exports the interpolated fields to the output directory defined in [`Export.ExportFolder`](parameters/Export/ExportFolder.md) as `MeanFlow.h5`.

### 2. Without interpolation

When `needInterpolation: false`, FELiCS assumes that the imported mean flow was previously exported by FELiCS (typically as `MeanFlow.h5`) on the same computational mesh as the current run.

Instead of interpolating the fields, FELiCS directly maps the imported degrees of freedom onto that of the current FEM spaces.

This avoids the interpolation step and can significantly reduce preprocessing time for large 3D simulations.

Example configuration:

```json
{
    "Case": {
        "needInterpolation": false,
        "MeshFilePath": "mesh.msh"
    },

    "FlowInput": {
        "MeanFlowFilePath": "MeanFlow.h5"
    }
}
```

## `.fel` file structure

A `.fel` file stores mean-flow information using the HDF5 format.

The file should contain:
- the coordinates of the points where the mean flow is defined,
- and the corresponding mean-flow variables.

Vector variables are stored component-wise by appending the coordinate name to the vector name.

Example structure for a 2D Cartesian flow containing coordinates `(x,y)`, velocity vector `u`, and density `rho`:

```text
/MeanFlow/x
/MeanFlow/y
/MeanFlow/ux
/MeanFlow/uy
/MeanFlow/rho
```

> 💡Note: the name of coordinates (and vector components) must be coherent with the coordinate system defined in [Case.CoordinateSystem](parameters/Case/CoordinateSystem.md). 
> E.g. `(x, y, z)` for 3D cartesian coordinates or `(x, r, t)` for 3D cylindrical coordinates.

## `.h5` file structure
The structure of this files is described in the FELiCS output files.

<!-- # Mean flow file
FELiCS processes the mean flow information with a .fel file using the .h5 format. It contains the mesh, mean flow field as well as the forcing and response domains.

Here is a python script example of how the mean flow file can be written: 
``` python
import  h5py

case = {
    'DATA_FILE':        'input_base_flow/data.mat', #input base flow data
    'MESH_FILE':        'input_mesh/data.mat', # input mesh
    'SAVE_FILE':        'myFELiCScase/meanflow.fel', #saving directory
}

# The data from the mesh and base flow are loaded

hf = h5py.File(case['SAVE_FILE'], 'w')
hf.create_dataset('/MeanFlow/x', data=x)
hf.create_dataset('/MeanFlow/r', data=y)
hf.create_dataset('/MeanFlow/ux', data=ux)
hf.create_dataset('/MeanFlow/uy', data=uy)
hf.create_dataset('/MeanFlow/uz', data=uz)
hf.create_dataset('/MeanFlow/nulam', data=nu_mol)

# One can include here an eddy viscosity field
hf.create_dataset('/MeanFlow/nuturb', data=nut)

# One can include here a sponge field
hf.create_dataset('/MeanFlow/spg', data=spg)
    
hf.create_dataset('/MeanFlow/responseDomain', data=Wresponse)
    hf.create_dataset('/MeanFlow/forcingDomain', data=Wforcing)
hf.close()
```

# Mean flow file

The base flow can be obtained from numerical simulations, experimental results or analytical models. A RANS mean field can also be computed using the finite element Newton solver `FlowSolver.py` integrateed in FELiCS.
The relevant mean flow information (velocities, pressure, viscosity, forcing and response domains, ...) are encapsulated in a **.fel** file. 
The construction of such a file is presented in [fel file in FELiCS](fel_file.md). -->