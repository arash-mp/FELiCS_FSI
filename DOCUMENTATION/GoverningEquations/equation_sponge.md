# Sponge functions

The sponge functions introduce additional source terms into the governing equations to stabilize numerical calculations. They can be applied to any variable—whether vector or scalar—within the system of equations being solved. In FELiCS, sponge functions are incorporated by iterating over all equations where sponge functionality is enabled and adding the corresponding damping terms to each.

Motivation:
- a sponge region can prevent nonphysical reflections and stabilizing the solution.
- a sponge region introduces additional damping terms to the linearized equations, ensuring that perturbations decay in the outflow or far-field, thus mimicking open boundary conditions and improving spectral accuracy.

Use case:
- in some open-domain cases, to supress nonphysical reflections.
- in cases with artificial boundaries to damp spurious oscillations or recirculation zones near the domain edges.
- in linear stability analysis on cases with non-reflecting boundary conditions to decay eigenmodes or response fields smoothly

References:
- [Mani 2012](https://doi.org/10.1016/j.jcp.2011.10.017)
- [Colombo et al. 2016](https://doi.org/10.1016/j.compfluid.2016.09.019)

Nomenclature:
- $\sigma$: sponge coefficient, defined by user
- $\phi$: state variable, can be any vector or scalar field that is solved
- $\phi_{target}$: The target distribution for state variable, needs to be defined by user
- $X_\phi$: test function for the equation regarding the state varable

Assumptions:
- fluctuations at sponge region are ignored
- reflections at sponge region are damped

## Nonlinear equations
When sponge is activated, the following term is added to the left hand side of the equation regarding to the state variable.
$$
-\sigma (\phi -\phi_{target})
$$

### Weak form
The weak form of the nonlinear sponge term, as implemented in FELiCS is
$$
\int_\Omega -j \sigma (\overline{\phi} - \phi_{target})\cdot X_{\phi} dx
$$

## Linear equations
The linear form of the sponge term is
$$
-\sigma \phi'
$$

### Weak form
The weak form of the linearized sponge term as implemented in FELiCS, is
$$
\int_\Omega -j\sigma \phi ' \cdot X_\phi dx
$$