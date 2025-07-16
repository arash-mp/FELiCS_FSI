# Sponge functions

The sponge functions introduce additional source terms into the governing equations to stabilize numerical calculations. They are mainly used as an additional boundary treatment for outflow or far-field boundaries. In FELiCS, if a sponge coefficient is given, sponge functions are incorporated by iterating over all equations and adding the corresponding damping terms to each.

Motivation:

A sponge region introduces additional damping terms to the linearized equations, ensuring that perturbations decay in the outflow or far-field, thus mimicking open boundary conditions and improving spectral accuracy.

Where to be used:
- In far-field and outflow regions, to supress nonphysical reflections.
- In cases with artificial boundaries to damp spurious oscillations or recirculation zones near the domain edges.
- In linear stability and resolvent analysis: As an "extended" boundary condition, combined with a zero-Drichlet condition, to decay eigenmodes or response fields smoothly.

References:
- [Mani 2012](https://doi.org/10.1016/j.jcp.2011.10.017)
- [Colombo et al. 2016](https://doi.org/10.1016/j.compfluid.2016.09.019)

Nomenclature:
- $\sigma$: sponge coefficient, defined by user. Should be zero in the main parts of the domain where a physical solution is desired, and be greater than zero in the respective boundary region.
- $\phi_{i}$: state variables for which the equations are solved (e.g.: fluctuations/modes)
- $\phi_{i, target}$: target value for the state variables in the boundary region; assumed to be zero for linear problems;  needs to be defined by the user if a non-linear problem is solved.
- $X_{\phi_{i}}$: FEM test function for the corresponding equation 

Assumptions:
- The fluctuations in sponge regions are damped. 

## Nonlinear equations
When the sponge is activated, the following term is added to the nonlinear equations for all state variables $\phi_{i}$:
$$
-\sigma (\phi_{i} -\phi_{i, target})
$$

### Weak form
The weak form of the nonlinear sponge term, as implemented in FELiCS is
$$
\int_\Omega -j \sigma (\overline{\phi_{i}} - \phi_{i,target})\cdot X_{\phi_{i}} dx
$$

## Linear equations
The linear form of the sponge term is
$$
-\sigma \phi_{i}'
$$

### Weak form
The weak form of the linearized sponge term, as implemented in FELiCS, is
$$
\int_\Omega -j\sigma \phi_{i} ' \cdot X_{\phi_{i}} dx
$$
