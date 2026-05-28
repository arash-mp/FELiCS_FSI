# FELiCS configuration file

FELiCS runs are controlled by a main configuration file, written in the `json` format.

This file defines both the parameters of the run (i.e. numerical and physical options) as well as the paths to the other FELiCS [input files](../index.md/#input-files).

The configuration file is directly passed to FELiCS from the command line with the `-f` (or `--file`) option:

```bash
FELiCS -f <configuration_file.json>
```

> 💡 **Note:** The configuration file is expected to be located in the running directory of the case, and the paths to other FELiCS files are interpreted relative that directory.

## File format

The configuration file uses the [JSON](https://www.json.org/json-en.html) format.

Some important JSON rules are:
- strings must use double quotes,
- boolean values are written as `true` or `false`,
- comments are not allowed,
- trailing commas are not allowed.

For example:
```JSON
{
    "Case": {
        "AnalysisMode":         "Resolvent",
        "CoordinateSystem":     "Cartesian",
        "nDim":                 2,
        "needInterpolation":    true
    }
}
```

## FELiCS parameters

The `settings.json` file is organized into 6 **parameter categories**. 

Each category groups parameters related to one aspect of the simulation setup. Use the following links to access the parameters definition of each category:

| Category | Purpose |
| - | - |
| [BoundaryCondition](parameters/BoundaryCondition/index.md) | *Path to the boundary-condition file* |
| [Case](parameters/Case/index.md) | *Physical model, geometry, analysis type, mesh, and equation setup* |
| [Export](parameters/Export/index.md) | *Output directory and export options* |
| [FlowInput](parameters/FlowInput/index.md) | *Mean-flow input settings* |
| [IOResolvent](parameters/IOResolvent/index.md) | *Resolvent and input-output analysis settings* |
| [Numerics](parameters/Numerics/index.md) | *Discretization and solver-related parameters* |

A typical configuration file will look like:
```
{
    "BoundaryCondition": {},
    "Case": {},
    "Export": {},
    "FlowInput": {},
    "IOResolvent": {},
    "Numerics": {}
}
```

## Defaults and parameter handling

FELiCS does not read arbitrary parameters from the configuration file. When a run starts, FELiCS initializes its internal dictionary of supported parameters and default values. The values provided in `<your_settings>.json` then overwrite these defaults.

This has two important consequences:
- parameters omitted from `<your_settings>.json` keep their default value;
- parameters that are not registered in FELiCS are ignored.

The list of supported parameters and their defaults is defined in the {py:func}`getAllSettingsDict() <FELiCS.parameters.config.getAllSettingsDict>` method of the FELiCS source code, in:

```bash
src/FELiCS/parameters/config.py
```

> 💡 **Note for developers:** When adding a new parameter to FELiCS, remember to add it to the default parameter dictionary in `config.py`; otherwise, it will not be read from `<your_settings>.json`.

## FELiCS-specific formatting conventions
Some FELiCS inputs uses additional conventions from JSON standards.

### Real/complex valued numbers

Real-valued parameters can be given as standard JSON numbers:
```JSON
"Omegas": [1.0, 2.0, 5.0e-1]
```

However, complex-valued parameters are given as strings:
```JSON
"EigenValueGuess": ["1.0+0.1j", "1.0-0.1j"]
```

Real and complex values can also be combined when the parameter accepts a list:
```JSON
"Omegas": [1.0, "2.0+0.5j"]
```

### Paths

File paths are interpreted relative to the directory from which FELiCS is executed, unless absolute paths are provided.

## Example

The following example is the configuration file from [tutorial 2: Modal Analysis](../../Tutorials/modal_analysis.md):

```json
{
    "BoundaryCondition": {
        "BCsFilePath":              "bc_modal.json"
    },
    "Case": {
        "AnalysisMode":             "Modal",
        "CoordinateSystem":         "Cartesian",
        "m":                        0.0,
        "MeshFilePath":             "cylinder_wake.msh",
        "MolVisc":                  0.02,
        "MolViscModel":             "Constant",
        "nDim":                     2,
        "Reaction":                 false,
        "SetOfEquations": {
            "Momentum": {
                "Equation":         "NSPrimitive",
                "Variable":         "u"
            },
            "Mass": {
                "Equation":         "Continuity",
                "Variable":         "p"
            },
            "Energy": {
                "Equation":         "None",
                "Variable":         "None"
            },
            "Species": {
                "Equation":         "None",
                "Variable":         "None"
            },
            "EquationOfState": {
                "Equation":         "None",
                "Variable":         "None"
            }
        },
        "SpeciesFilePath":          "",
        "TurbulenceModel":          "None"
    },
    "Export": {
        "ExportFolder":             "output_dir"
    },
    "FlowInput": {
        "MeanFlowFilePath":         "base_flow_for_FELiCS.fel"
    },
    "IOResolvent": {
        "ForcingBoundaryIndices":   [],
        "ForcingMode":              "Body",
        "Omegas":                   [],
        "ResponseCoeff":            []
    },
    "Numerics": {
        "EigenValueGuess":          ["0.7"],
        "nSolut":                   100,
        "PolynomialOrder": {
            "u":                    2,
            "p":                    1
        }
    }
}
```

<!-- 💡 **Note that FELiCS offers two approaches for importing the mean flow:**

1. **With interpolation** (default): In general, the meanflow is given on a different grid than that of the linear computation (defined in .msh mesh file). The meanflow variables stored in a .fel file are then (linearly) interpolated onto the computational grid. For large meshes, particularly 3D cases, this can be **time-consuming**. The interpolated meanflow is then exported to `MeanFlow.h5`. Corresponding inputs are:
    ```json
    {
        "needInterpolation": true,
        "MeshFilePath": "your_mesh.msh",
        "FlowInput": {"MeanFlowFilePath": "your_meanflow.fel"}
    }
    ```

2. **Without interpolation**: if you already have an interpolated meanflow corresponding to your computational mesh (e.g. from a previous run), you may skip the interpolation alltogether by indicating the FELiCS exported `MeanFlow.h5` as input. This is particularly useful when running multiple linearized analyses with the same meanflow and mesh. Corresponding inputs are:
    ```json
    {
        "needInterpolation": false,
        "MeshFilePath": "your_mesh.msh",
        "FlowInput": {"MeanFlowFilePath": "MeanFlow.h5"}
    }
    ```

## Specific formatting for FELiCS settings and json
* **booleans** are defined in lowercase (`true`, `false`)
* **`EigenValueGuess`** and **`Omegas`** can be either
  * list of strings for complex values: `["1.0-1j", "1.0+1j"]`
  * list of floats for real values: `[1, 2, 1.2e-1]`
  * a combination of both: `["1.0-1j", 1, "1.0+1j", 2]`
* json accepts exponential notation for all float-type inputs

>**Warning:** only the parameters listed by `config.getAllSettingsDict()` in the `src/FELiCS/parameters/config.py` file will be considered by FELiCS. Add your new parameters there to be able to use them in the code. -->
