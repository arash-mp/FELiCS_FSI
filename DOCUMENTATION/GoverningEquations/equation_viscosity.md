# Viscosity models
Depending on the flow state, the effective viscosity $\mu_\textrm{eff} = \mu_m + \mu_t$ has to be treated with different models. In laminar, incompressible flows, the turbulent eddy viscosity $\mu_t = 0$ and only the molecular viscosity $\mu_m$ remains. For incompressible, isothermal flows, the molecular viscosity can usually be assumed to be spatially constant. For compressible and/or flows with a variable temperature field, the molecular viscosity is often assumed to be a function of the temperature according to Sutherland's law. In turbulent flows, $\mu_t \neq 0$ and an eddy viscosity model such as the classical Boussinesq model is required. Both of these models are briefly described in the following.

## Sutherland model
The Sutherland model describe how molecular viscosity changes with temperature for a specific gas.

Motivation:
- The molecular viscosity of a specific gas can change with temperature, the Sutherland model can describe it.

Use case:
- cases with large gradient temperature which can affect molecular viscosity significantly. 

References:
- [Sutherland et al. 1893](https://doi.org/10.1080/14786449308620508)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $T$: Temperature
- $T_0$: Reference temperature
- $\mu$: molecular dynamic viscosity at temperature $T$
- $\mu_0$: Reference dynamic viscosity at reference temperature $T_0$
- $C$: Sutherland's constant

Assumptions/justification:
- Ideal gas behavior
- Negligible intermolecular forces except during collisions
- Viscosity depends only on temperature

The Sutherland model is
$$
\mu = \mu_0 \frac{T_0 +C}{T+C}(\frac{T}{T_0})^{3/2}
$$
It is used for calculating dynamic molecular viscosity of the mean flow.

__NOTE: The linearized form of the Sutherland model is not implemented in FELiCS yet.__


## Boussinesq model
The Boussinesq model is a turbulence closure approach that approximates the effect of turbulent fluctuations by introducing an eddy viscosity, which relates the Reynolds stresses to the mean strain rate. Eddy viscosity represents the enhanced momentum transport due to turbulence and is typically modeled as an additional scalar in diffusion term in momentum equation. 

Motivation:
- to closure the RANS equaitons, Reynolds stress meeds to be modeled 

Use case:
- most turbulent cases using RANS equations

References:
- [Boussinesq 1872](https://gallica.bnf.fr/ark:/12148/bpt6k56673076/f2.item.texteImage)

Nomenclature:
- $\rho \overline{\mathbf{u}'\mathbf{u}'}$: Reynolds stress tensor
- $\mu_t$: eddy viscosity
- $\mathbf{u}$: velocity vector
- $k$: turbulent kenetic energy
- $\mathbf{I}$: identity tensor


Assumptions/justification:
- Turbulent Reynolds stresses are proportional to the mean strain rate.
- Turbulent mixing is analogous to molecular viscosity, using an effective eddy viscosity.
- Turbulent viscosity is isotropic.

The Boussinesq assumption models the Reynolds stress like 
$$
-\rho \overline{\mathbf{u}'\mathbf{u}'} =\mu_t (\nabla +\nabla^T )\mathbf{u} - \frac{2}{3} k \mathbf{I}
$$
combine its first term on right hand side with diffusion term in Momentum equation, the viscosity can be expressed like:
$$
\mu_{Eff} = \mu + \mu_t
$$

### Frozen eddy viscosity
In this method, the fluctuation of eddy visicosity is ignore, which means
$$
\mu' = 0
$$

### Linear equations
There are various methods on linearization eddy viscosity. The following methods are implemented in present version of FELiCS:
- Frozen eddy viscosity

