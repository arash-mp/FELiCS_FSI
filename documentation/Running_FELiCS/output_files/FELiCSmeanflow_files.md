# Mean-flow output files

FELiCS exports interpolated mean-flow fields as:

```text
MeanFlow.h5
```

The file uses the HDF5 format and stores the mean-flow quantities used internally by FELiCS.

The exported fields are defined on the FELiCS [export mesh](FELiCSmesh_files.md). 

The file is exported in the output directory defined by the [Export.ExportFolder](../input_files/parameters/Export/ExportFolder.md) parameter.


## Contents

The exported quantities are not restricted to the FELiCS state vector.

The file may contain:
- flow variables,
- thermodynamic quantities,
- turbulent viscosity fields,
- sponge regions,
- forcing and response domains,
- and additional fields required by the selected physical models.

## Variable naming convention

The naming convention follows the same structure as the `.fel` input files.

Scalar variables are stored directly using their variable name.

Vector variables are stored component-wise by appending the coordinate name to the vector name.

Example:

```text
/x
/y

/ux
/uy

/rho
```

See also:
- [Mean-flow input files](../input_files/mean_files.md)