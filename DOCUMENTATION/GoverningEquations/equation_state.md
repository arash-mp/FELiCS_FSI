# Equations of state
The equations of state are required when thermodynamic processes become relevant in the flow field. The equations of state relate a set of state variables to each other.

References:
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $p$: pressure
- $p^{(0)}$: slow pressure component in low-Mach expansion
- $T$: temperature
- $\rho$: density
- $R$: specific gas constant

## Ideal gas equation
Example use case:
- Compressible flow including acoustics: [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Assumptions:
- Valid for a hypothetical *ideal* gas.

The ideal gas equation is
$$
p = \rho R T
$$

### Linearized ideal gas equation
The linearized form of the ideal gas equation is
$$
p' = R (\rho' \overline{T} +\overline{\rho} T')
$$


## Low-Mach equation
Assumptions:
- Ideal gas with constant mean pressure, i.e. mean pressure in the entire flow field is constant due to low Mach numbers and the pressure field is set to the slow component of the pressure $p^{(0)}$ fixed at some point of the domain, *e.g.* at the outlet.
- Density only changes with temperature, but not with pressure.

Example use case:
- Reacting turbulent jet flame: [Input-output analysis](./../Tutorials/input_output_analysis.md)

The low-Mach equation is
$$
\rho = \frac{p^{(0)}}{RT}
$$

### Linearized low-Mach equation
The linearized form of the low-Mach equation is
$$
\rho' =  -\frac{\overline{\rho}T'}{\overline{T}}
$$