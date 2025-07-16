# Equations of state
Motivation:
- equations for relating state variables to each other
- important to close equation system

Use case:
- flows with acoustics, flows with combustion/chemical reactions
- tutorial example: turbulent flame? heat conduction?

References:
- [Demange2024 AIAA](https://doi.org/10.2514/6.2024-3170)
- [Kaiser2023 CnF](https://doi.org/10.1016/j.combustflame.2023.112778)

Nomenclature:
- $p$: pressure
- $p_0$: inlet pressure
- $T$: temperature
- $\rho$: density
- $R$: specific gas constant

## Ideal gas
Assumptions/justification:
- ideal gas assumption: for a specific gas, its density, pressure and temperature can effect each other with the gas constant
### Nonlinear equations
$$
p = R \rho T
$$
### Linear equations
The linear form of the ideal gas equation is:
$$
p' = R (\rho' \overline{T} +\overline{\rho} T')
$$
## Low-Mach (ideal gas with constant pressure)
Assumptions/justification:
- mean pressure in the entire flow field is constant due to low Mach numbers and the pressure field is set to the inlet pressure $p_0$

### Nonlinear equations
$$
\rho = \frac{p_0}{RT}
$$
### Linear equations
The linear form of the Low-Mach ideal gas equaiton is:
$$
\rho' =  -\frac{\overline{\rho}T'}{\overline{T}}
$$
