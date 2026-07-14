# Case parameters

The `Case` category contains the main parameters defining the physical problem, mathematical formulation, and type of linear analysis performed by FELiCS.

## Contents

| Parameter | Description |
| - | - |
| [`AnalysisMode`](AnalysisMode.md) | Type of linear analysis performed by FELiCS |
| [`CalculateAdjoint`](CalculateAdjoint.md) | Enables the computation of adjoint modes in modal analysis |
| [`CoordinateSystem`](CoordinateSystem.md) | Coordinate system used for the computation |
| [`needInterpolation`](needInterpolation.md) | Controls how mean-flow fields are imported onto the FELiCS computational mesh |
| [`MeshFilePath`](MeshFilePath.md) | Path to the computational mesh file |
| [`MixtureFilePath`](MixtureFilePath.md) | Path to the mixture-definition file |
| [`MolVisc`](MolVisc.md) | Constant molecular viscosity value |
| [`MolViscModel`](MolViscModel.md) | Molecular viscosity model used in the simulation |
| [`MolViscPerturbModel`](MolViscPerturbModel.md) | Model used for fluctuating molecular viscosity |
| [`Reaction`](Reaction.md) | Enables reaction-related quantities in FELiCS |
| [`SpeciesFilePath`](SpeciesFilePath.md) | Path to the species-definition file (currently unused, reserved for upcoming developments) |
| [`m`](m.md) | Wavenumber in the spectral direction |
| [`nDim`](nDim.md) | Number of spatial dimensions resolved by the mesh |
| [`PrandtlNumber`](PrandtlNumber.md) | Global Prandtl number used for thermal diffusion |
| [`SetOfEquations`](SetOfEquations.md) | Defines the transported equations and variables |
| [`TurbulenceModel`](TurbulenceModel.md) | Turbulence model used in the simulation |

## Example
```json
"Case": {
    "AnalysisMode": "Resolvent",
    "CoordinateSystem": "Cartesian",
    "nDim": 2,
    "MeshFilePath": "mesh.msh"
}
```
## Related documentation

- See the [mesh file](../../gmsh_files.md), [mixture file](../../mix_files.md), and [theory pages](../../../../GoverningEquations/index.md) for more background on the physical and numerical models configured through this category.

## Detailed content
```{toctree}
:maxdepth: 1

AnalysisMode.md
CalculateAdjoint.md
CoordinateSystem.md
needInterpolation.md
MeshFilePath.md
MixtureFilePath.md
MolVisc.md
MolViscModel.md
MolViscPerturbModel.md
Reaction.md
SpeciesFilePath.md
m.md
nDim.md
PrandtlNumber.md
SetOfEquations.md
TurbulenceModel.md

```