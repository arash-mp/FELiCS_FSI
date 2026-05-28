# FELiCS mixture files

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