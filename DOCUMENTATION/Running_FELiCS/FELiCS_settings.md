# Configuration files
The input parameters are loaded from the three following **.json** files:
- `settings.json`: contains most of the information about the FELiCS run.
- `boundaries.json`: contains specific information about the boundary conditions.
- `mixture.json`(optional): contains specific information about the physical properties of the fluid considered.

## Structure of the `settings.json` file
The typical `settings.json` file is divided into 6 main sections:


``` json
{
"BoundaryCondition":{
    "BCsFilePath":      (str)   path to the "boundaries.json" file,
},
"Case":{
    "AnalysisMode":     (str)   type of analysis, one of: ("Modal", "Resolvent", "InputOutput"),
    "CalculateAdjoint": (bool)  if we compute the adjoint spectrum in "Modal" analysis,
    "CoordinateSystem": (str)   coordinate system, one of: ("Cylindrical", "Cartesian"),
    "m":                (float) wavenumber of fluctuations in the third, spectral dimension (for a 2D grid),
    "MeshFilePath":     (str)   path to the mesh file.
    "MixtureFilePath":  (str)   path to the "mixture.json" file.
    "MolVisc":          (float) kinematic viscosity (dynamic viscosity in compressible flows) value when "MolViscModel" is set to "Constant"
    "MolViscModel":     (str)   type of laminar viscosity model used for the viscosity,
    "MolViscPerturbModel": (dict)   type of laminar viscosity model used for fluctuating laminar viscosity, specifying [type, empirical constants], e.g.:
        {
            "type":"Sutherland mean",
            "Constants": {"Ts":170.672}
        },
    "nDim":             (int)   number of dimensions resolved in the mesh.
    "PrandtlNumber":    (int)   spatially global Prandtl number, only relevant for cases where thermal diffusion or conduction is involved
    "Reaction":         (bool)  if using chemistry reactions (DEPRECATED? => DEFINED IN THE MIXTURE FILE?)
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
    "ForcingNorm":      (str)   norm type for the forcing term ("Resolvent" analysis), one of: ("TKE" - default, "Chu")
    "ResponseNorm":     (str)   norm type for the response term ("Resolvent" analysis), one of ("TKE" - default, "Chu)
    "Omegas":           (list)  angular frequency (see formatting in previous section),
},
"Numerics":{
    "EigenValueGuess":  (list)  eigenvalue guesses for "Modal" analysis (see formatting in previous section),
    "nSolut":           (int)   number of solutions to compute (number of nearest eigenvalues in  "Modal" analysis, number of resolvent modes per frequency in "Resolvent" analysis)
    "PolynomialOrder":  (dict)  sets the polynomial order for each transported variables (list must match "SetOfEquations"). e.g.: {"u": 2,"T": 1,"rho": 1}
}}
```
>**Note:** The name of the file is not important. 
The `settings.json` is directly given to FELiCS when running from the command line via the `-file` flag:
```bash
FELiCS -file settings.json
```

## Structure of the `boundaries.json` file
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
The name of this json file should be set in the main settings file via the variable "BCsFilePath".


## Structure of the `mixture.json` file
>**TODO: ** review (Thomas?)

In most cases that do not involve chemistry modelling or species transport, the `mixture.json` is not needed. 

Below is the structure of the `mixture.json` file for a case using chemistry:

```json
{
"Species":                        (dict) list of species that are considered, providing how they are calculated and additional species properties, e.g.:
    { 
        "phi":{                   (dict) name of species
            "calc":"transported", (str) specify species that are returned, "all": all species, "transported": only species with separate transport equations, "constraint": only passive species without transport equations,
            "Sc":0.9              (float) Schmidt number of the species
        }
    },
"Reaction_mechanism": (dict) list with type of reaction model and all required settings, e.g.:
{
    "type": "KaiserCnF2023", (str) type of reaction model (currently only "KaiserCnF2023" is implemented)
    "additional_fields":["prefactor"], (str) additional fields to be read in from the input *.fel file 
    "reactions":[{ (dict) list of reactions that provides all involved educts and products and their corresponding stoichiometric coefficients
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
