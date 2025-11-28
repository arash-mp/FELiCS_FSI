# Species transport equation
The species transport equation describes the transport behavior of a species k of a mixture.

Assumptions:
- source terms need to be defined based on a given model, e.g. [reaction model](./equation_reaction.md)
- Schmidt number, which relates mass diffusivity with eddy viscosity, is often assumed to be constant

Example use case:
- passive scalar transport: [Kaiser et al. 2021](https://doi.org/10.1017/jfm.2021.151)
- reacting turbulent jet flame: [Input-output analysis](./../Tutorials/input_output_analysis.md)

References:
- [Kaiser et al. 2021](https://doi.org/10.1017/jfm.2021.151)
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Nomenclature:
- $Y_k$: mass fraction of species k
- $\mathbf{u}$: velocity vector
- $\rho$: density
- $\breve{D}_\textrm{eff}$: density-premultiplied mass diffusivity, effective
- $\breve{D}$: density-premultiplied mass diffusivity, molecular
- $\breve{D}_t$: density-premultiplied mass diffusivity, turbulent
- $f_{Y_k}$: source term (e.g. chemical reaction rate due to a flame)
- $\mathrm{Sc}$: Schmidt number
- $X_{Y_k}$: test function for species transport equation



### Nonlinear species transport equation
The nonlinear species transport equation is
$$
\rho \frac{\partial Y_k}{\partial t} + \rho \mathbf{u} \cdot \nabla Y_k = \nabla \cdot (\breve{D}_\textrm{eff} \nabla Y_k) + f_{Y_k}
$$
where $\breve{D}_\textrm{eff}$ is the effective density-premultiplied mass diffusivity, and $f_Y$ is a source term, which occurs for example due to a chemical reaction of a flame. The effective density-premultiplied mass diffusivity consists of the molecular and turbulent diffusivity defined as $\breve{D}_\textrm{eff} = \breve{D} + \breve{D}_t$. In laminar flows $\breve{D}_t = 0$. In turbulent flows, one approach is to link the turbulent diffusivity $\breve{D}_t$ to the turbulent eddy viscosity $\mu_t$ via $\breve{D}_t = \mu_t / \mathrm{Sc}_t$, where $\mathrm{Sc}_t = \mu_t / \breve{D}_t$ is the turbulent Schmidt number that is often assumed to be constant. For more details on turbulent eddy viscosity, see [Turbulence models](./equation_turbulence_models.md).

```{note}
The diffusive term is correct only for very simple models, such as passive scalar transport and reaction progress variables. To generalize this, the diffusion should be built with the gradient of the molar fraction $X_k$ (See eq. (1.45) in Theoretical and Numerical Combustion, Poinsot & Veynante)
```

```{note}
The nonlinear equations are currently *not* implemented in FELiCS.
```

```{note}
The density-premultiplied effective mass diffusivity is defined as $\breve{D}_\textrm{eff} = \rho D_\textrm{eff}$, with $[D_\textrm{eff}] = \mathrm{m}^2/\mathrm{s}$ if using SI units.
```



### Linearized species transport equation
The linear species transport equation is
$$
\overline{\rho} \frac{\partial Y_k'}{\partial t} + \rho ' \overline{\mathbf{u}} \cdot \nabla \overline{Y}_k + \overline{\rho} \mathbf{u}' \cdot \nabla \overline{Y}_k + \overline{\rho} \, \overline{\mathbf{u}} \cdot \nabla Y_k' = \nabla \cdot \left(\breve{D}_\textrm{eff}' \nabla \overline{Y}_k\right) + \nabla \cdot \left(\overline{\breve{D}}_\textrm{eff} \nabla Y_k'\right) + f_{Y_k}'
$$


#### Weak form of the linearized species transport equation
The weak form of the linearized species transport equation with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega \omega \overline{\rho} \hat{Y}_k X_{Y_k}^* \mathrm{d}\mathbf{x} = \int_\Omega \mathrm{j} \, \overline{Y}_k \, \nabla \cdot \left( \hat{\rho} \, \overline{\mathbf{u}} \, X_{Y_k}^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} \mathrm{j} \, \overline{Y}_k \, \hat{\rho} \, (\overline{\mathbf{u}} \cdot \mathbf{n}) X_{Y_k}^* \, \mathrm{d}\mathbf{s} + \int_\Omega \mathrm{j} \, \overline{Y}_k \, \nabla \cdot \left( \overline{\rho} \, \hat{\mathbf{u}} \, X_{Y_k}^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} \mathrm{j} \, \overline{Y}_k \, \overline{\rho} \, (\hat{\mathbf{u}} \cdot \mathbf{n}) X_{Y_k}^* \mathrm{d}\mathbf{s} + \int_\Omega \mathrm{j} \, \hat{Y}_k \, \nabla \cdot \left( \overline{\rho} \, \overline{\mathbf{u}} \, X_{Y_k}^* \right) \, \mathrm{d}\mathbf{x} - \int_{\partial \Omega} \mathrm{j} \, \hat{Y}_k \, \overline{\rho} \, (\overline{\mathbf{u}} \cdot \mathbf{n}) X_{Y_k}^* \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} \, \hat{\breve{D}}_\textrm{eff} \, \nabla \overline{Y}_k \cdot \nabla X_{Y_k}^* \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} \, \hat{\breve{D}}_\textrm{eff} \, \left(\nabla \overline{Y}_k \cdot \mathbf{n}\right) X_{Y_k}^* \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} \, \overline{\breve{D}}_\textrm{eff} \, \nabla \hat{Y}_k \cdot \nabla X_{Y_k}^* \, \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} \, \overline{\breve{D}}_\textrm{eff} \, \left(\nabla \hat{Y}_k \cdot \mathbf{n}\right) X_{Y_k}^* \, \mathrm{d}\mathbf{s} + \int_\Omega f_{Y_k}' \, X_{Y_k}^* \, \mathrm{d}\mathbf{x}
$$

```{note}
The density-premultiplied mass diffusivity is defined as $\breve{D}_\textrm{eff} = \rho D_\textrm{eff}$, with $[D_\textrm{eff}] = \mathrm{m}^2/\mathrm{s}$ if using SI units.
```
