# Reaction models
docs


## Eddy break up model
docs

Motivation:

docs

Where to be used:
- points

References:
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Nomenclature:
- $\dot{\Omega}$: Source term in energy equation results from chemical reaction
- $A_{EBU}$: Model constant
- $\varepsilon$: dissipation ratio of turbulent kinetic energy
- $k$: turbulent kinetic energy
- $\rho$: density
- $c$: progress variable, which is a scalar field


Assumptions:
- docs

## Nonlinear equations
The source term in energy equation results from the chemical reaction is modeled as
$$
\dot{\Omega} = A_{EBU} \rho c (1-c)
$$

## Linear equations
The linear form of the eddy break up model equation is:
$$
\dot{\Omega}' = A_{EBU} (\rho' (\overline{c} - \overline{c}^2) + \overline{\rho}(c'-2\overline{c}c'))
$$
where the fluctions of $k$ and $\varepsilon$ are ignected.
