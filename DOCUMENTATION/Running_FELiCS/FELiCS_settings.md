# Setting files
The input parameters are loaded from the three following **.json** files:
- `settings.json`: contains most of the information about the FELiCS run.
- `boundaries.json`: contains specific information about the boundary conditions.
- `mixture.json`: contains specific information about the physical properties of the fluid considered.

## Structure of the `settings.json` file
The typical `settings.json` file is divided into 6 main sections:

``` json
{
"BoundaryCondition":{
    "BCsFilePath":      (str)   path to the "boundaries.json" file,
},
"Case":{
    "AnalysisMode":     (str)   type of analysis ("Modal", "Resolvent", "InputOutput"),
    "CalculateAdjoint": (bool)  if we compute the adjoint spectrum in "Modal" analysis,
    "CoordinateSystem": (str)   coordinate system ("Cylindrical", "Cartesian"),
    "m":                (float) wavenumber of fluctuations in the Fourier dimension (for 2D analysis),
    "MeshFilePath":     (str)   path to the mesh file.
    "MixtureFilePath":  (str)   path to the "mixture.json" file.
    "MolVisc":          (float) laminar viscosity (dynamic viscosity in compressible!) value when "MolViscModel" is set to "Constant"
    "MolViscModel":     (str)   type of laminar viscosity model used for laminar viscosity,
    "nDim":             (int)   number of dimensions resolved in the mesh.
    "Reaction":         (bool)  if using chemistry reactions (DEPRECATED?)
    "SetOfEquations":{  (dict)  with fields: [equations, variants, transported variables], all defined as str. e.g.:  
        "Momentum":{
            "Equation":"NSPrimitive",
            "Variable":"u"
            }
        },
    "SpeciesFilePath":  (str)   path to the SpeciesFile (DEPRECATED?),
    "TransVelFluc":     (bool)  if fluctuations are defined in the Fourier dimension (DUPLICATED with "m"?)
    "TurbulenceModel":  (str)   type of turbulence model (currently "None", "Constant", or "File")
},
"Export":{
    "ExportFolder":     (str)   directory (relative to running dir) where outputs are saved,
    "Video":            (bool)  if xdmf files used to generate a video of the mode are exported, 
},
"FlowInput":{
    "AveragingDirection":(str)  dimension in which the mean flow is averaged, set to "None" if unused (DEPRECATED?)
    "MeanFlowFilePath": (str)   mean flow file path and name to load.
},
"IOResolvent":{
    "ForcingBoundaryIndices": (list of int) indices of boundaries where forcing is applied ("InputOutput" analysis),
    "ForcingCoeff":     (list of int)   variables onto which forcing is applied ("InputOutput" analysis),
    "ForcingMode":      (str)   type of forcing "Body" or "Boundary" ("InputOutput" analysis)
    "ForcingNorm":      (str)   norm type for the forcing term ("Resolvent" analysis)
    "ResponseNorm":     (str)   norm type for the response term ("Resolvent" analysis)
    "Omegas":           (list)  angular frequency (see formatting in previous section),
},
"Numerics":{
    "EigenValueGuess":  (list)  eigenvalue guesses for "Modal" analysis (see formatting in previous section),
    "nSolut":           (int)   number of solutions to compute (in "Modal" and "Resolvent" analysis)
    "PolynomialOrder":  (dict)  sets the polynomial order for each transported variables (list must match "SetOfEquations"). e.g.: {"u": 2,"T": 1,"rho": 1}
}}
```
>**Note:** The name of the file is not important. 
The `settings.json` is directly given to FELiCS when running from the command line via the `-file` flag:
```bash
FELiCS -file settings.json
```

## Structure of the `boundaries.json` file
The structure of this file is expected to evolve soon. The current structure is:
```json
{
"var1" (str):
    [{"ID": (int), "type": (str) "Dirichlet" or "Neumann", "value": (float)},
"var2" (str):
    [{"ID": (int), "type": (str) "Dirichlet" or "Neumann", "value": (float)},
}
```
where:
* `"var" (str)` is the name of the variable considered. 
* `"ID": (int)` is the index of the boundary considered. For gmsh *.msh* meshes, this must correspond to the index of an existing *PhysicalNames* entry.
* `"type": (str)` is the type of BC applied at the boundary for the variable considered. Currently only accepts `"Dirichlet"` or `"Neumann"`.
* `"value": (float)` is the value imposed on the variable (or its gradient). **Always 0.0 ??**

A typical `boundaries.json` file is:

```json
{
"300":
    [{"variable": "ux", "type": "Dirichlet", "value": 0.0},
    {"variable": "uy", "type": "Dirichlet", "value": 0.0},
    {"variable": "p", "type": "Neumann", "value": 0.0}],
"301":
    [{"variable": "ux", "type": "Neumann", "value": 0.0},
    {"variable": "uy", "type": "Neumann", "value": 0.0},
    {"variable": "p", "type": "Dirichlet", "value": 0.0}]
}
```

## Structure of the `mixture.json` file
>**Warning:** No idea how this is organized. Please, Thomas or someone who knows about this completes the documentation here.

In most cases that do not involve chemistry modelling to describe the fluid, the `mixture.json` is not used and a dummy file is passed instead. The structure of the dummy file is:

```json
{
"species":{
    "Air" :{
        "calc":"constraint"
    }
},
"reaction":{}
}
```

The following is an example of the `mixture.json` file for a case using chemistry:

```json
{
"Species":{
    "progress":{
        "calc":"transported",
        "Sc":0.9
    }
},
"Pr":0.9,
"Viscosity":{
    "type":"Constant",
    "Constants":{
        "Viscosity":1.0
    }
},
"Reaction_mechanism":{
    "type": "KaiserCnF2023",
    "additional_fields":["prefactor"],
    "reactions":[{
        "educts": [],
        "stochiometricCoefficientsEducts": [],
        "products": ["progress"],
        "stochiometricCoefficientsProducts": [1.0]
    }]
}}
```

## Specific formatting for FELiCS settings and json
* **booleans** are defined in lowercase (`true`, `false`)
* **`EigenValueGuess`** and **`Omegas`** can be either
  * list of strings for complex values: `["1.0-1j", "1.0+1j"]`
  * list of floats for real values: `[1, 2, 1.2e-1]`
  * a combination of both: `["1.0-1j", 1, "1.0+1j", 2]`
* json accepts exponential notation for all float-type inputs

>**Warning:** only the parameters listed by `config.getAllSettingsDict()` in the `src/FELiCS/parameters/config.py` file will be considered by FELiCS. Add your new parameters there to be able to use them in the code.
