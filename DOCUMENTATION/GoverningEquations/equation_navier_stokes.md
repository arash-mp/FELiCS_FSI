# Navier--Stokes equations

The Navier–Stokes equations represent the fundamental principles of fluid mechanics, expressing the conservation of momentum and mass in a fluid. They provide a comprehensive description of fluid flow behavior across a wide range of physical situations. Because these equations capture the essential physics of Newtonian fluids--such as air, water, and gases under 'standard' conditions--they form the fundamental basis for analyzing and predicting flow dynamics.

Example use case:
- incompressible: [Solving for the base flow of the cylinder wake](./../Tutorials/cylinder_wake.md)

References:
- [Barkley et al. 2006](https://doi.org/10.1209/epl/i2006-10168-7)
- [Towne et al. 2018](https://doi.org/10.1017/jfm.2018.675)
- [Müller et al. 2024](https://doi.org/10.1017/jfm.2024.679)
- [Demange et al. 2024](https://doi.org/10.2514/6.2024-3170)

Nomenclature:
- $\mathbf{u}$: velocity vector
- $p$: pressure
- $\rho$: density
- $\mathbf{\tau}$: viscous stress tensor
- $\mu$: dynamic viscosity
- $\mathbf{I}$: identity tensor
- $X_p$: test function for mass equation
- $\mathbf{X}_\mathbf{u}$: test function for momentum equations


## Incompressible
Assumptions:
- primitive variables
- $\rho = \textrm{const}$, this is reasonable if Mach number is low, no large temperature gradients, no acoustics, no combustion
- no additional source terms
- no gravitational forces

### Nonlinear Navier--Stokes equations
The mass equation is
$$
\nabla \cdot \mathbf{u} = 0
$$

The momentum equations are
$$
\rho \frac{\partial \mathbf{u}}{\partial t} + \rho(\mathbf{u}\cdot \nabla) \mathbf{u} + \nabla p - \nabla \cdot \mathbf{\tau} = 0
$$

where the viscous stress tensor $\mathbf{\tau}$ is
$$
\mathbf{\tau} = \mu_\textrm{eff}(\nabla + \nabla ^T)\mathbf{u}
$$
and where $\mu_\textrm{eff} = \mu + \mu_t$ is the effective dynamic viscosity (see [Viscosity models](./equation_viscosity.md) for details).

### Mean flow equations
We consider the flow field to be comprised of a  time-invariant base flow, which can be either a time-averaged flow or fixed point solution (base flow), and the perturbation, such that 

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\mathbf{q}'(\mathbf{x},t), \qquad \mathbf{q}=(\mathbf u,p)^{T}
$$

Inserting this into the Navier--Stokes equations and taking the time-average we get the base  flow equations.

The mass equation is
$$
\nabla \cdot \overline{\mathbf{u}} = 0
$$

The momentum equations are
$$
\rho(\overline{\mathbf{u}}\cdot \nabla) \overline{\mathbf{u}} + \nabla \overline{p} - \nabla \cdot \overline{\mathbf{\tau}} = 0
$$

where the mean viscous stress tensor $\overline{\mathbf{\tau}}$ is
$$
\overline{\mathbf{\tau}} = \overline{\mu}_\textrm{eff}(\nabla + \nabla ^T)\overline{\mathbf{u}}
$$

#### Weak form of the mean flow equations
```{note}
The weak form of the incompressible nonlinear base flow equations is only implicitly included in the generalized compressible form of the Navier--Stokes equations (see below) and is, therefore, not stated here.
```





### Linearized Navier--Stokes equations
The mass equation is
$$
\nabla \cdot \mathbf{u}' = 0
$$

The momentum equations are
$$
\rho \frac{\partial \mathbf{u}'}{\partial t} + \rho [(\overline{\mathbf{u}}\cdot \nabla) \mathbf{u}' + (\mathbf{u}' \cdot \nabla) \overline{\mathbf{u}}] + \nabla p' - \nabla \cdot \mathbf{\tau}' = 0
$$

where the fluctuating viscous stress tensor $\tau'$ is
$$
\mathbf{\tau}' = \mu_\textrm{eff}'[(\nabla + \nabla ^T)\overline{\mathbf{u}}]+ \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T) \mathbf{u}']
$$
and where $\mu_\textrm{eff}'$ is the fluctuating effective dynamic viscosity (see [Viscosity models](./equation_viscosity.md) for details).

#### Weak form of the linearized Navier--Stokes equations
```{note}
The weak form of the incompressible linearized Navier--Stokes equations is only implicitly included in the generalized compressible form of the Navier--Stokes equations (see below) and is, therefore, not stated here.
```









## Compressible
Assumptions:
- primitive variables
- no additional source terms
- no gravitational forces

### Nonlinear Navier--Stokes equations
The mass equation is
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
$$

The momentum equations are
$$
\rho \frac{\partial \mathbf{u}}{\partial t} + \rho(\mathbf{u}\cdot \nabla) \mathbf{u} + \nabla p - \nabla \cdot \tau = 0
$$

where the viscous stress tensor $\tau$ is
$$
\tau = \mu_\textrm{eff}[(\nabla + \nabla ^T)\mathbf{u} - \frac{2}{3} (\nabla \cdot \mathbf{u})\mathbf{I}]
$$
where $\mu_\textrm{eff}$ is the effective dynamic viscosity (see [Viscosity models](./equation_viscosity.md) for details) and where $\mathbf{I}$ is the identity tensor.

### Mean flow equations

We consider the flow field to be comprised of a  time-invariant base flow, which can be either a time-averaged flow or fixed point solution (base flow), and the perturbation, such that 

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\mathbf{q}'(\mathbf{x},t), \qquad \mathbf{q}=(\mathbf{u},p,\rho)^{T}
$$

```{note}
For the compressible equtions, the overbar represents a favre-abverage for the velocity and a Reynolds average for pressure and density.  
```

Inserting this into the Navier--Stokes equations and taking the average we get the base  flow equations

The mass equation is
$$
\nabla \cdot (\overline{\rho}\, \overline{\mathbf{u}}) = 0
$$

The momentum equations are
$$
\overline{\rho}(\overline{\mathbf{u}}\cdot \nabla) \overline{\mathbf{u}} + \nabla \overline{p} - \nabla \cdot \overline{\tau} = 0
$$


where the mean viscous stress tensor $\overline{\mathbf{\tau}}$ is
$$
\overline{\tau} = \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T)\overline{\mathbf{u}} - \frac{2}{3} (\nabla \cdot \overline{\mathbf{u}})\mathbf{I}]
$$





### Weak form of the mean flow equations
The weak form of the nonlinear mean flow equations, as implemented in FELiCS, is
$$
\int_\Omega \mathrm{j} \nabla \cdot \left(\mathbf{X}_\mathbf{u}^* \otimes \overline{\rho}\overline{\mathbf{u}} \right) \cdot \overline{\mathbf{u}} \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} \mathrm{j} \overline{\rho}\left( \left(\overline{\mathbf{u}} \otimes \mathbf{X}_\mathbf{u}^*\right) \cdot \overline{\mathbf{u}} \right)\cdot \mathbf{n} \, \mathrm{d}\mathbf{s}  + \int_\Omega \mathrm{j} \overline{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} \mathrm{j} \overline{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} \overline{\mathbf{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} \overline{\mathbf{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s} = 0
$$
where the mean viscous stress tensor $\overline{\tau}$ is
$$
\mathbf{\overline{\tau}} = \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T)\overline{\mathbf{u}}]
$$
and where $\overline{\mu}_\textrm{eff}$ is the mean effective dynamic viscosity (see [Viscosity models](./equation_viscosity.md) for details).


### Linearized Navier--Stokes equations
The mass equation is
$$
\frac{\partial \rho'}{\partial t} + \nabla \cdot (\overline{\rho} \mathbf{u}' + \rho' \overline{\mathbf{u}}) = 0
$$

The momentum equations are
$$
\overline{\rho} \frac{\partial \mathbf{u}'}{\partial t} + \overline{\rho} [(\overline{\mathbf{u}}\cdot \nabla) \mathbf{u}' + (\mathbf{u}' \cdot \nabla) \overline{\mathbf{u}}] + \rho' (\overline{\mathbf{u}}\cdot \nabla)\overline{\mathbf{u}} + \nabla p' - \nabla \cdot \tau ' = 0
$$

where the fluctuating viscous stress tensor $\tau'$ is

$$
\tau' = \mu_\textrm{eff}'[(\nabla + \nabla ^T)\overline{\mathbf{u}} - \frac{2}{3} \nabla \cdot \overline{\mathbf{u}} \mathbf{I}]+ \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T) \mathbf{u}' - \frac{2}{3}( \nabla \cdot \mathbf{u}' )\mathbf{I}]
$$
and where $\mu_\textrm{eff}'$ is the effective dynamic viscosity (see [Viscosity models](./equation_viscosity.md) for details).

#### Weak form of the linearized Navier--Stokes equations
The weak form of the linearized Navier--Stokes equations with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega \omega \hat{\rho}X_p^* \mathrm{d}\mathbf{x} = \int_\Omega \mathrm{j} \widehat{\rho\mathbf{u}} \cdot \nabla X_p^* \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} \mathrm{j} \widehat{\rho\mathbf{u}} \cdot \mathbf{n} X_p^* \, \mathrm{d}\mathbf{s}
$$
and
$$
\int_\Omega \omega \overline{\rho} \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} = - \int_\Omega \mathrm{j} \left(\overline{\rho}\overline{\mathbf{u}} \cdot \nabla \right) \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} -\int_\Omega \mathrm{j} \left(\overline{\rho}\hat{\mathbf{u}} \cdot \nabla \right) \overline{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} - \int_\Omega \mathrm{j} \left(\hat{\rho}\overline{\mathbf{u}} \cdot \nabla \right) \overline{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} + \int_\Omega \mathrm{j} \hat{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} \mathrm{j} \hat{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega \mathrm{j} \hat{\mathbf{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} \mathrm{j} \hat{\mathbf{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s}
$$
with
$$
\widehat{\rho\mathbf{u}} = \overline{\rho}\mathbf{\hat{u}} + \hat{\rho}\mathbf{\overline{u}}
$$
```{note}
In the linearized equations, the convective term is *not* integrated by parts.
```
