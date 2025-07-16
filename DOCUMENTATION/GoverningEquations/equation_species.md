# Species transport equation
The species transport equation describes the transport behavior of a species $Y$ of a mixture, or the transport of a passive scalar.

Motivation:
- Component concentration can affect physical properties, so it needs to be track when a mixture is applied;
- The transport of passive scalar fields are very important in some cases. For example progress variables in combustion models.

Use case:
- in flows with chemical reactions
- in flows with mixture which has non-uniform concentration distribution
- tutorial example: turbulent flame?

References:
- [Kaiser et al. 2021](https://doi.org/10.1017/jfm.2021.151)
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Nomenclature:
- $Y$: species or passive scalar
- $\mathbf{u}$: velocity vector
- $\rho$: density
- $\breve{D}_\textrm{eff}$: density-premultiplied mass diffusivity, effective
- $\breve{D}$: density-premultiplied mass diffusivity, molecular
- $\breve{D}_t$: density-premultiplied mass diffusivity, turbulent
- $f_Y$: source term (e.g. chemical reaction rate due to a flame)
- $\mathrm{Sc}$: Schmidt number
- $X_Y$: test function for species transport equation

Assumptions:
- The density $\rho$ usually changes due to chemical reactions, changes of component concentrations, high-Mach numbers, etc...
- Source terms can be defined based on reaction / combustion models
- often Schmidt number is assumed to be constant
- no gravitational forces

### Nonlinear species transport equation
The nonlinear species transport equation is
$$
\rho \frac{\partial Y}{\partial t} + \rho \mathbf{u} \cdot \nabla Y = \nabla \cdot (\breve{D}_\textrm{eff} \nabla Y) + f_{Y}
$$
where $\breve{D}_\textrm{eff}$ is the effective density-premultiplied mass diffusivity, and $f_Y$ is a source term, which occurs for example due to a chemical reaction of a flame. The effective density-premultiplied mass diffusivity consists of the molecular and turbulent diffusivity defined as $\breve{D}_\textrm{eff} = \breve{D} + \breve{D}_t$. In laminar flows $\breve{D}_t = 0$. In turbulent flows, one approach is to link the turbulent diffusivity $\breve{D}_t$ to the turbulent eddy viscosity $\mu_t$ via $\breve{D}_t = \mu_t / \mathrm{Sc}_t$, where $\mathrm{Sc}_t = \mu_t / \breve{D}_t$ is the turbulent Schmidt number that is often assumed to be constant. For more details on turbulent eddy viscosity, see [Viscosity models](./equation_viscosity.md).

```{note}
The nonlinear equations are currently *not* implemented in FELiCS.
```

```{note}
The density-premultiplied effective mass diffusivity is defined as $\breve{D}_\textrm{eff} = \rho D_\textrm{eff}$, with $[D_\textrm{eff}] = \mathrm{m}^2/\mathrm{s}$ if using SI units.
```



### Linearized species transport equation
The linear species transport equation is
$$
\overline{\rho} \frac{\partial Y'}{\partial t} + \rho ' \overline{\mathbf{u}} \cdot \nabla \overline{Y} + \overline{\rho} \mathbf{u}' \cdot \nabla \overline{Y} + \overline{\rho} \overline{\mathbf{u}} \cdot \nabla Y' = \nabla \left(\cdot \breve{D}_\textrm{eff}' \nabla \overline{Y}\right) + \nabla \cdot \left(\overline{\breve{D}}_\textrm{eff} \nabla Y'\right) + f_{Y}'
$$


#### Weak form of the linearized species transport equation
The weak form of the linearized species transport equation with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega \omega \overline{\rho} \hat{Y} X_Y^* \mathrm{d}\mathbf{x} = \int_\Omega j \, \overline{Y} \, \nabla \cdot \left( \hat{\rho} \, \overline{\mathbf{u}} \, X_Y^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} j \, \overline{Y} \, \hat{\rho} \, (\overline{\mathbf{u}} \cdot \mathbf{n}) X_Y^* \, \mathrm{d}\mathbf{s} + \int_\Omega j \, \overline{Y} \, \nabla \cdot \left( \overline{\rho} \, \hat{\mathbf{u}} \, X_Y^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} j \, \overline{Y} \, \overline{\rho} \, (\hat{\mathbf{u}} \cdot \mathbf{n}) X_Y^* \mathrm{d}\mathbf{s} + \int_\Omega j \, \hat{Y} \, \nabla \cdot \left( \overline{\rho} \, \overline{\mathbf{u}} \, X_Y^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} j \, \hat{Y} \, \overline{\rho} \, (\overline{\mathbf{u}} \cdot \mathbf{n}) X_Y^* \mathrm{d}\mathbf{s} - \int_\Omega j \, \hat{\breve{D}}_\textrm{eff} \, \nabla \overline{Y} \cdot \nabla X_Y^* \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \, \hat{\breve{D}}_\textrm{eff} \, \left(\nabla \overline{Y} \cdot \mathbf{n}\right) X_Y^* \, \mathrm{d}\mathbf{s} - \int_\Omega j \, \overline{\breve{D}}_\textrm{eff} \, \nabla \hat{Y} \cdot \nabla X_Y^* \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \, \overline{\breve{D}}_\textrm{eff} \, \left(\nabla \hat{Y} \cdot \mathbf{n}\right) X_Y^* \, \mathrm{d}\mathbf{s} + \int_\Omega f_Y' \, X_Y^* \, \mathrm{d}\mathbf{x}
$$

```{note}
The density-premultiplied mass diffusivity is defined as $\breve{D}_\textrm{eff} = \rho D_\textrm{eff}$, with $[D_\textrm{eff}] = \mathrm{m}^2/\mathrm{s}$ if using SI units.
```