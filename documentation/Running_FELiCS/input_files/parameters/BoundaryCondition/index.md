# BoundaryCondition parameters

The `BoundaryCondition` category contains parameters related to the boundary-condition definition used by FELiCS.

*For now, this category only includes the path to the boundary condition file used in the run.*

## Contents

| Parameter | Description |
| - | - |
| [`BCsFilePath`](BCsFilePath.md) | Path to the boundary-condition file |

## Example

```json
{
    "BoundaryCondition": {
        "BCsFilePath": "boundaries.json"
    }
}
```

## Related documentation

- [Boundary-condition file](../../bc_files.md)

## Detailed content
```{toctree}
:maxdepth: 1

BCsFilePath.md

```

