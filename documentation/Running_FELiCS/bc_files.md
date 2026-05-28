# Boundary-condition file

The boundary-condition file defines the boundary conditions applied by FELiCS on the computational domain boundaries.

The file is provided in JSON format and referenced in the main configuration file through [`BoundaryCondition.BCsFilePath`](../Configuration/BoundaryCondition/BCsFilePath.md).

## Boundary identifiers

Each entry of the boundary-condition file corresponds to one mesh boundary.

The dictionary keys (`"1"`, `"2"`, ...) should match the *Physical Group numbers* defined in the mesh file.

For example, the following `.msh` section:

```text
$PhysicalNames
6
1 1 "Inlet"
1 2 "Symmetry"
1 3 "Outlet"
1 4 "Top"
1 5 "Wall"
2 6 "all"
$EndPhysicalNames
```

defines:
- boundary `1` → `"Inlet"`
- boundary `2` → `"Symmetry"`
- boundary `3` → `"Outlet"`
- ...

In the boundary-condition file, these identifiers are referenced through their numerical indices:

```json
{
    "1": {...},
    "2": {...},
    "3": {...}
}
```

The Physical Group names (`"Inlet"`, `"Wall"`, ...) are currently ignored by FELiCS.

## Structure

Each boundary entry contains:
- a boundary-condition name,
- and optionally a list of boundary-condition specifications.

General structure:

```json
{
    "1": {
        "name": "wall"
    }
}
```

## Supported boundary-condition types

| Name | Description |
| - | - |
| `"zeroDirichlet"` | Imposes homogeneous Dirichlet conditions on all variables |
| `"wall"` | Wall boundary condition |
| `"custom"` | User-defined boundary conditions |
| `"symmetry"` | Symmetry boundary condition |
| `"none"` | No boundary condition applied |

## Custom boundary conditions

The `"custom"` and `"symmetry"` boundary types require a `"specifics"` entry defining the imposed conditions.

Each specification contains:
- `"variable"`: variable name,
- `"type"`: boundary-condition type,
- `"value"`: imposed value.

| Type | Description |
| - | - |
| `"Dirichlet"` | Imposes the variable value |
| `"Neumann"` | Imposes the variable gradient |
| `"None"` | No condition applied |

## Example

```json
{
    "1": {
        "name": "zeroDirichlet"
    },

    "2": {
        "name": "symmetry",
        "specifics": [
            {
                "variable": "ux",
                "type": "Dirichlet",
                "value": 0.0
            },
            {
                "variable": "uy",
                "type": "Neumann",
                "value": 0.0
            },
            {
                "variable": "p",
                "type": "Dirichlet",
                "value": 0.0
            }
        ]
    },

    "3": {
        "name": "wall"
    }
}
```

<!-- # FELiCS boundary condition files

The structure of this file is:
```json
{
"ID1":{            (int as a str)    ID of the boundary. For gmsh *.msh* meshes, this must correspond to the index of an existing *PhysicalNames* entry.
    "name":        (str)             name of boundary type; one of: ("zeroDirichlet", "wall", "custom", "symmetry","none")
    "specifics":{  (dict)            only for boundary types "custom" and "symmetry"
        "variable":(str)             the name of the variable considered; has to be one of the state vector varaibles, defined in the "SetOfEquations" part in the general config.json file
        "type"    :(str)             type of boundary conditions, one of: ("Dirichlet", "Neumann", "None")
        "value"   :(float)           value imposed on the variable (or its gradient)
}}}
```

Here is an expamle for a `boundaries.json` file:

```json
{
    "1": {
        "name": "zeroDirichlet"
    },
    "2": {
        "name": "symmetry",
        "specifics": [
            {
                "variable": "ux",
                "type": "Dirichlet",
                "value": 0.0
            },
            {
                "variable": "uy",
                "type": "Neumann",
                "value": 0.0
            },
            {
                "variable": "p",
                "type": "Dirichlet",
                "value": 0.0
            }
        ]
    },
    "3": {
        "name": "wall"
    }
}
```
>**Note:** The name of the file is not important. 
The name of this json file should be set in the main settings file via the variable "BCsFilePath". -->