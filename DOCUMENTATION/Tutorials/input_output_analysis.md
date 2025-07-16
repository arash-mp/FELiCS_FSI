# Tutorial 4: Input/Output analysis

## Goal of the Tutorial
The goal of this tutorial is to provide a step-by-step guide on performing input/output analysis in FELICS for a reacting flow. By the end of this tutorial, you will be able to:

- Configure reaction mechanisms and species transport in FELICS.
- Write customed boundary conditions.
- Perform input/output analysis with boundary forcing for reactive flows using the mass, momentum, species, and energy equations.

## Case Definition
In this tutorial, we will perform a compressible input/output analysis about a 2D reacting jet flame base flow. The analysis focuses on understanding how the system responds to a given boundary forcing. The axial component of the velocity $u_x$ is harmonicaly forced at the inlet. The case uses a progress variable approach to model the reaction and includes species transport effects.

The 3D configuration is axisymmetric, allowing us to solve the 3D problem for a given azimuthal wavenumber and work with a 2D mesh. The reaction is modeled using a progress variable that tracks the reaction progress from reactants to products.

### Mesh generation
We generate the 2D mesh on GMSH for the combustor geometry. The mesh file `KITBurnerWallSep.msh` is used in this tutorial, which represents a confined burner configuration.

**Warning:** For 2D computations, FELiCS handles only triangular elements only.

The mesh should include proper boundary identification for:
- Inlet boundaries (for fuel/air injection)
- Outlet boundaries (for product extraction)
- Wall boundaries (no-slip conditions)
- Symmetry axis (for axisymmetric cases)

### Base Flow
This tutorial case folder is located in ```felics2.0/TUTORIALS/TURB_FLAME```.

Our base flow is a time-averaged reacting flow field that includes:
- Velocity components ($u_x, u_r, u_\theta$)
- Pressure field $p$
- Species concentration (progress variable $c$)
- Temperature field $T$ (derived from progress variable)
- Turbulent viscosity field $\nu_t$

The base flow is stored in the file ```KIT_confined.fel```. This flow field represents the steady-state solution of the reacting flow equations and serves as the base state around which we perform the input/output analysis. The base flow axial velocity is displayed in [Figure1](#UXMean) 
![Figure1](../TUTORIALS/input_ouput_tutorial/pic/MeanFlow.png) <a id="fig:UXMean"></a>

Figure 1: Mean flow axial velocity

## Input/Output Analysis Parameters

### Boundary conditions
Here we set the axisymmetric boundary conditions in the [```boundaries.json```](../../TUTORIALS/TURB_FLAME/boundaries.json) file.

The boundary conditions for reacting flows include additional considerations for species transport (here the progress $c'$):

| Boundary | $u'_x$ | $u'_r$ | $u'_\theta$ | $p'$ | $c'$|
|:----------|:-----------|:-----------|:-----------|:-----------|:-----------|
| <code style="color : Darkorange">Forcing (Boundary 1)</code> | None | Dirichlet | Dirichlet | None | Dirichlet |
| <code style="color : Darkorange">Outlet (Boundary 2)</code> | None | None | None | Dirichlet | None |
| <code style="color : Darkorange">Symmetry (Boundary 3)</code> | Neumann | Dirichlet | Neumann | Neumann | Neumann |
| <code style="color : Darkorange">Walls (Boundaries 4,5,6)</code> | Dirichlet | Dirichlet | Dirichlet | Neumann | Neumann |

**Note:** The forcing boundary (Boundary 1) is where external perturbations are applied to study the system's response. The progress variable boundary conditions ensure proper species transport at each boundary.

### Reaction Mechanism
The reaction mechanism is defined in the [```Mixture.json```](../../TUTORIALS/TURB_FLAME/Mixture.json) file. For this tutorial, we use a simplified progress variable approach:

```json
{
"Species":{
    "progress":{
        "calc":"transported",
        "Sc":0.9
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

This configuration defines:
- A transported progress variable with Schmidt number 0.9
- A reaction mechanism based on the KaiserCnF2023 model
- Additional fields for reaction prefactor calculations
- A single reaction that produces the progress variable

### Settings
The setting file [```turb_flame.json```](../../TUTORIALS/TURB_FLAME/turb_flame.json) contains all the analysis information. The explanation of each field is provided in [Setting files](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/development/DOCUMENTATION/Running_FELiCS/FELiCS_settings.md?ref_type=heads).

Key settings for input/output analysis with reacting flows:

#### Coordinate System and Geometry
```json
"CoordinateSystem": "Cylindrical",
"m": 0.0,
"nDim": 2
```
The coordinates system is set to Cylindrical with azimuthal wavenumber m=0 (axisymmetric mode).

#### Analysis Mode
```json
"AnalysisMode": "Input-Output",
"CalculateAdjoint": true
```
This enables input/output analysis with adjoint calculation for sensitivity analysis.

#### Equations Set
```json
"SetOfEquations": {
    "Momentum": {
        "Equation": "NSPrimitive",
        "Variable": "u"
    },
    "Mass": {
        "Equation": "Continuity",
        "Variable": "p"
    },
    "Species": {
        "Equation": "Non-conservative",
        "Variable": "progress"
    },
    "Energy": {
        "Equation": "ProgressVariableLinear",
        "Variable": "None"
    },
    "EquationOfState": {
        "Equation": "Low-Mach",
        "Variable": "None"
    }
}
```
This configuration solves:
- Momentum equations for velocity fluctuations
- Mass conservation for pressure
- Species transport for progress variable
- Energy equation coupled to progress variable
- Low-Mach equation of state allows to assume hte mean ressure to be constant. The density is only impacted by the temperature.

#### Physical Properties
```json
"MolVisc": 0.0002,
"MolViscModel": "Constant",
"PrandtlNumber": 0.9,
"TurbulenceModel": "File"
```
- Molecular viscosity set to 0.0002
- Prandtl number for thermal diffusion
- Turbulent viscosity read from file

#### Forcing Configuration
```json
"IOResolvent": {
    "ForcingBoundaryIndices": [1],
    "ForcingCoeff": [0, 1],
    "ForcingMode": "Boundary",
    "Omegas": [314]
}
```
This configures:
- Forcing applied at boundary index 1
- Forcing coefficients for different variables
- Boundary forcing mode
- Analysis frequency (ω = 314)

#### Reaction Settings
```json
"Reaction": false,
"MixtureFilePath": "Mixture.json"
```
While direct reaction is disabled, the mixture file provides the reaction mechanism for the progress variable approach.

## Running the analysis
Your case folder should look like:
```bash
.
├── mesh/
│   └── KITBurnerWallSep.msh
├── meanFlow/
│   └── KIT_confined.fel
├── turb_flame.json
├── boundaries.json
├── Mixture.json
└── Out/
```

Run the analysis with the command:
```bash
FELiCS -f turb_flame.json
```

The input/output analysis will:
1. Load the reacting flame base flow
2. Set up the linearized equations for momentum, mass, species, and energy
3. Apply forcing at the specified boundary
4. Solve for the system response and compute transfer functions
5. Calculate adjoint fields for sensitivity analysis
6. Export results to the output directory

The analysis provides insights into:
- How the reacting flow responds to external perturbations
- The coupling between fluid dynamics and chemical reactions
- Sensitivity of the flame to different forcing mechanisms
- Transfer functions between input forcing and output response


