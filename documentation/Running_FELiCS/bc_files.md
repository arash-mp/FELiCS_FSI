# FELiCS boundary condition files

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