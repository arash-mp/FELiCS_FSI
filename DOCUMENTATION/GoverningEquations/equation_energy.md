# Energy equations
The energy equations are usually required for compressible flows when thermodynamic changes (e.g. changes in enthalpy or temperature) in the flow field become relevant, such as flows featuring acoustics or combustion/chemical reactions. Various forms of the energy equation exist, each tailored to specific physical problems and the assumptions made during their derivation. The motivation behind using energy equations lies in their role as fundamental expressions of energy conservation (or its equivalent forms) in fluid flows involving thermodynamic processes.

References:
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $\mathbf{u}$: velocity vector
- $p$: pressure
- $T$: temperature
- $T_b$: burnt, adiabatic flame temperature
- $T_u$: unburnt temperature
- $\rho$: density
- $\gamma$: heat capacity ratio
- $\kappa$: thermal conductivity
- $\mathbf{\tau}$: viscous stress tensor
- $\mu$: dynamic viscosity
- $c$: Progress variable
- $\mathbf{I}$: identity tensor
- $X_h$: test function for energy equation



## Energy pressure equation
Motivation?

Assumptions:
- primitive variables
- constant Prandtl number
- constant heat capacity ratio
- molecular viscosity adheres to Sutherland's law
- fluctuations in thermal conductivity are tied to the given Prandtl number and, thus, eddy viscosity

Example use case:
- Compressible flow including acoustics: [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

### Nonlinear energy pressure equation
The nonlinear energy pressure equation is
$$
\frac{\partial p}{\partial t} + \mathbf{u} \cdot \nabla p + \gamma p(\nabla \cdot \mathbf{u}) - (\gamma - 1)\left[ \nabla (\kappa \nabla T) + \tau : \nabla \mathbf{u}\right] = 0
$$

with the viscous stress tensor $\tau$ being defined as

$$
\tau = \mu \left[\left(\nabla + \nabla ^T\right)\mathbf{u} - \frac{2}{3} (\nabla \cdot \mathbf{u}) \mathbf{I}\right].
$$

```{note}
The nonlinear energy pressure equation is currently not implemented in FELiCS.
```

### Linearized energy pressure equation
The linear energy pressure equation is
$$
\frac{\partial p'}{\partial t} + \overline{\mathbf{u}} \cdot \nabla p' + \mathbf{u}' \cdot \nabla \overline{p} + \overline{\gamma}\left[ \overline{p}(\nabla \cdot \mathbf{u}') + p'(\nabla \cdot \overline{\mathbf{u}}) \right] - (\overline{\gamma} - 1)\left[ \nabla (\overline{\kappa} \nabla T' + \kappa' \nabla \overline{T} + \overline{\tau} : \nabla \mathbf{u}' + \tau' : \nabla \overline{\mathbf{u}}) \right] = 0,
$$

where the mean viscous stress tensor $\overline{\tau}$ is defined as

$$
\mathbf{\overline{\tau}} = \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T)\overline{\mathbf{u}} - \frac{2}{3} (\nabla \cdot \overline{\mathbf{u}}) \mathbf{I}].
$$

The fluctuating viscous shear stress tensor $\tau'$ is defined as

$$
\tau' = \mu' \left[\left(\nabla + \nabla ^T\right)\overline{\mathbf{u}} - \frac{2}{3} (\nabla \cdot \overline{\mathbf{u}}) \mathbf{I}\right] + \overline{\mu}\left[(\nabla + \nabla ^T) \mathbf{u}' - \frac{2}{3} (\nabla \cdot \mathbf{u}') \mathbf{I}\right]
$$

where $\overline{\mu}_\textrm{eff}$ and $\mu_\textrm{eff}'$ are the mean and fluctuating effective dynamic viscosity, respectively (see [Viscosity models](./equation_viscosity.md) for details).


#### Weak form of the linearized energy pressure equation
The weak form of the linearized energy pressure equation with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega \omega \hat{p} \, X_h^* \, \mathrm{d}\mathbf{x} = \int_\Omega \mathrm{j} \, \nabla \cdot \left( X_h^* \, \overline{\mathbf{u}} \right) \, \hat{p} \, \mathrm{d}\mathbf{x} + \int_\Omega \mathrm{j} \, \nabla \cdot \left( X_h^* \, \hat{\mathbf{u}} \right) \, \overline{p} \, \mathrm{d}\mathbf{x} + \int_\Omega \mathrm{j} \, \overline{\gamma} \, \nabla \left( \overline{p} \, X_h^* \right) \cdot \hat{\mathbf{u}} \, \mathrm{d}\mathbf{x} + \int_\Omega \mathrm{j} \, \overline{\gamma} \, \nabla \left( \hat{p} \, X_h^* \right) \cdot \overline{\mathbf{u}} \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} \mathrm{j} (\overline{\gamma}+1) \, \overline{p} \, \hat{\mathbf{u}} \cdot \mathbf{n} X_h^* \, \mathrm{d}\mathbf{s} - \int_{\partial\Omega} \mathrm{j} (\overline{\gamma}+1) \, \hat{p} \, \overline{\mathbf{u}} \cdot \mathbf{n} X_h^* \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} (\overline{\gamma} - 1) \, \nabla X_h^* \cdot \left( \overline{\kappa} \nabla  \hat{T}  \right) \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} (\overline{\gamma} - 1) \, \overline{\kappa} \left( \nabla \hat{T} \cdot \mathbf{n} \right) X_h^* \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} (\overline{\gamma} - 1) \, \nabla X_h^* \cdot \left( \hat{\kappa} \nabla  \overline{T}  \right) \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} (\overline{\gamma} - 1) \, \kappa' \left( \nabla \overline{T} \cdot \mathbf{n} \right) X_h^* \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} (\overline{\gamma} - 1) \, \left(\nabla \cdot \overline{\mathbf{\tau}}\right) \cdot \left( \hat{\mathbf{u}}  X_h^* \right) \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} (\overline{\gamma} - 1) \, \left(\overline{\mathbf{\tau}} \hat{\mathbf{u}}\right) \cdot \mathbf{n} X_h^* \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} (\overline{\gamma} - 1) \, \left( \nabla \cdot \hat{\mathbf{\tau}} \right) \cdot \left( \overline{\mathbf{u}}  X_h^* \right) \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} (\overline{\gamma} - 1) \, \left(\hat{\mathbf{\tau}} \overline{\mathbf{u}}\right) \cdot \mathbf{n} X_h^* \, \mathrm{d}\mathbf{s}
$$





## Progress variable equation
This equation is used to describe an active flame model, which incorporates fluctuations of the reaction rate to capture turbulent flame dynamics more accurately ([Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)). This model relates the progress variable directly to the reaction rate, enabling a simplified yet effective representation of active flames under turbulent conditions.

Assumptions:
- progress variable and density are directly related through algebraic equations, which obviates the use of an explicit energy equation $\rightarrow$ the species equation also acts as an energy conservation equation
- incompressible flow due to (1) low Mach number, (2) flame being acoustically compact
- adiabatic walls (no heat losses), adiabatic flame temperature is reached in the entire domain
- constant heat capacity and constant specific gas constants
- all assumptions result in temperature only being a function of the progress variable, $T = T_u + (T_b - T_u)c$
- [low Mach equation of state](./equation_state.md), pressure is constant

Example use case:
- Reacting turbulent jet flame: [Input-output analysis](./../Tutorials/input_output_analysis.md)


### Nonlinear progress variable equation
See [Species transport equation](./equation_species.md).


### Linearized progress variable equation
See [Species transport equation](./equation_species.md).