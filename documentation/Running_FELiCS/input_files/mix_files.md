# Mixture file

The mixture file defines species properties and reaction models used by FELiCS for reacting-flow simulations.

The file is provided in JSON format and referenced in the main configuration file through [`Case.MixtureFilePath`](parameters/Case/MixtureFilePath.md).

In simulations that do not involve:
- chemistry modelling,
- species transport,
- or reacting flows,

the mixture file is typically not required.

## Structure

The mixture file may contain:
- a list of species definitions,
- and an optional reaction mechanism definition.

## Species

The `"Species"` entry defines the species considered in the simulation together with their associated properties.

Typical entries include:
- `"calc"`: how the species is treated by FELiCS,
- `"Sc"`: Schmidt number of the species.

### Example

```json
"Species": {
    "phi": {
        "calc": "transported",
        "Sc": 0.9
    }
}
```

| `calc` value | Description |
| - | - |
| `"all"` | Includes all species |
| `"transported"` | Includes species with dedicated transport equations |
| `"constraint"` | Includes passive constrained species |

## Reaction mechanism

The `"Reaction_mechanism"` entry defines the chemistry model used by FELiCS.

### Example

```json
"Reaction_mechanism": {
    "type": "KaiserCnF2023",

    "additional_fields": [
        "prefactor"
    ],

    "reactions": [
        {
            "educts": [],
            "stochiometricCoefficientsEducts": [],

            "products": [
                "progress"
            ],

            "stochiometricCoefficientsProducts": [
                1.0
            ]
        }
    ]
}
```

| Parameter | Description |
| - | - |
| `"type"` | Reaction model used by FELiCS |
| `"additional_fields"` | Additional fields imported from the mean-flow file |
| `"reactions"` | List of chemical reactions |

Currently implemented reaction models:

| Type | Description |
| - | - |
| `"KaiserCnF2023"` | Kaiser combustion model |

## Complete example

```json
{
    "Species": {
        "phi": {
            "calc": "transported",
            "Sc": 0.9
        }
    },

    "Reaction_mechanism": {
        "type": "KaiserCnF2023",

        "additional_fields": [
            "prefactor"
        ],

        "reactions": [
            {
                "educts": [],
                "stochiometricCoefficientsEducts": [],

                "products": [
                    "progress"
                ],

                "stochiometricCoefficientsProducts": [
                    1.0
                ]
            }
        ]
    }
}
```

<!-- # FELiCS mixture files

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
``` -->