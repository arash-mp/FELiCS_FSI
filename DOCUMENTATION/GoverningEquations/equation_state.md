# Equations of state
Motivation:
- equations for relating state variables to each other
- important to close equation system

Use case:
- flows with acoustics, flows with combustion/chemical reactions
- tutorial example: turbulent flame? heat conduction?

References:
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $p$: pressure
- $p_0$: inlet pressure
- $T$: temperature
- $\rho$: density
- $R$: specific gas constant

## Ideal gas equation
Assumptions/justification:
- classical ideal gas assumption, for a specific gas, its density, pressure and temperature are related to each other via the gas constant

The ideal gas equation is
$$
p = R \rho T
$$

### Linearized ideal gas equation
The linearized form of the ideal gas equation is
$$
p' = R (\rho' \overline{T} +\overline{\rho} T')
$$


## Low-Mach equation
Assumptions/justification:
- ideal gas with constant mean pressure, i.e. mean pressure in the entire flow field is constant due to low Mach numbers and the pressure field is set to the inlet pressure $p_0$
- therefore, the density only changes with temperature, but not with pressure

$$
\rho = \frac{p_0}{RT}
$$

### Linearized low-Mach equation
The linear form of the low-Mach ideal gas equation is
$$
\rho' =  -\frac{\overline{\rho}T'}{\overline{T}}
$$