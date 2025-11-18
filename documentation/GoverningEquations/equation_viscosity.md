# Viscosity models

This section describes the current models of **molecular viscosity** implemented in FELiCS.

## Sutherland model for molecular viscosity
The Sutherland model describes how molecular viscosity changes with temperature for a specific gas.

Assumptions:
- ideal gas behavior
- negligible intermolecular forces except during collisions
- molecular viscosity is only a function of the gas composition and the temperature

Example use case:
- Compressible flow including acoustics: [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

References:
- [Sutherland et al. 1893](https://doi.org/10.1080/14786449308620508)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $T$: Temperature
- $T_0$: Reference temperature
- $\mu$: molecular dynamic viscosity at temperature $T$
- $\mu_0$: Reference dynamic viscosity at reference temperature $T_0$
- $C$: Sutherland's constant


The Sutherland model is
$$
\mu = \mu_0 \frac{T_0 +C}{T+C}(\frac{T}{T_0})^{3/2}
$$
It is used for calculating dynamic molecular viscosity of the mean flow.

```{note}
The linearized form of the Sutherland model is not implemented in FELiCS yet.
```
