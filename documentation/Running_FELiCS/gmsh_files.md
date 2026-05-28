# Mesh file

FELiCS imports computational meshes using the [Gmsh](https://gmsh.info/) `.msh` format.

Both 2D and 3D meshes are supported.

## Supported mesh format

Meshes should be exported from Gmsh (or compatible meshing tools) using:

- **msh Version 2**
- **ASCII format**

Other `.msh` versions are currently not guaranteed to be supported by FELiCS.

## Physical groups

FELiCS relies on Gmsh *Physical Groups* to identify the different boundaries and regions of the computational domain.

It is therefore essential to define the appropriate Physical Groups during mesh generation.

These groups are used internally by FELiCS to:
- apply boundary conditions,
- identify forcing boundaries.

FELiCS identifies boundaries using the **Physical Group numbers**, not their names.

For example:

```geo
Physical Curve(1) = {5};
Physical Curve(2) = {6};
```

In this case:
- boundary `1` corresponds to curve `5`,
- boundary `2` corresponds to curve `6`.

The numerical identifiers (`1`, `2`, ...) are the values used in:
- the [`IOResolvent.ForcingBoundaryIndices`](../Configuration/IOResolvent/ForcingBoundaryIndices.md) parameter,
- the [boundary-condition](bc_files.md) file.

## Example

An example of mesh generation and usage in FELiCS is provided in [Tutorial 2](../Tutorials/modal_analysis.md).