# Navier--Stokes equations

The Navier--Stokes equations comprise the momentum and mass/continuity equations.

Motivation:
- fundamental equations for EVERY type flow
- conservation of momentum and mass

Use case:
- in every type of flow comprising Newton fluids (like air, water, gas flows at `normal' conditions)
- tutorial example incompressible with link: cylinder wake
- tutorial example compressible with link: cylinder wake

References:
- [Mueller et al. 2024](https://doi.org/10.1017/jfm.2024.679)
- [Barkley et al. 2006](https://doi.org/10.1209/epl/i2006-10168-7)
- [Towne et al. 2018](https://doi.org/10.1017/jfm.2018.675)
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
- $\rho = \textrm{const}$, this is reasonable if Mach number is low, no large temperature gradients/changes, no acoustics, no combustion
- no additional source terms
- no gravitational forces

### Nonlinear equations
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
and where $\mu_\textrm{eff}$ is the effective dynamic viscosity (see viscosity models for details).

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




#### Weak form
To solver the base flow equations with FELiCS we need to write it in weak form, reading
$$
\int_\Omega j \mathbf{\overline{u}} \cdot \nabla X_p^* \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \mathbf{\overline{u}} \cdot \mathbf{n} X_p^* \, \mathrm{d}\mathbf{s} = 0
$$

and

$$
\int_\Omega j \nabla \cdot \left(\mathbf{X}_\mathbf{u}^* \otimes \mathbf{u} \right) \cdot \mathbf{u} \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \left( \left(\mathbf{\overline{u}} \otimes \mathbf{X}_\mathbf{u}^*\right) \cdot \mathbf{\overline{u}} \right)\cdot \mathbf{n} \, \mathrm{d}\mathbf{s} + \int_\Omega j \overline{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \overline{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega j \mathbf{\overline{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \mathbf{\overline{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s} = 0
$$
where the mean viscous stress tensor $\overline{\tau}$ is
$$
\mathbf{\overline{\tau}} = \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T)\overline{\mathbf{u}}]
$$
and where $\overline{\mu}_\textrm{eff}$ is the mean effective dynamic viscosity (see viscosity models for details).




### Linear equations
The continuity equation is
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
and where $\mu_\textrm{eff}'$ is the fluctuating effective dynamic viscosity (see viscosity models for details).

#### Weak form
The weak form of the linearized Navier--Stokes equations with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega j \hat{\mathbf{u}} \cdot \nabla X_p^* \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \hat{\mathbf{u}} \cdot \mathbf{n} X_p^* \, \mathrm{d}\mathbf{s} = 0
$$
and
$$
\int_\Omega \omega \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} = - \int_\Omega j \left(\overline{\mathbf{u}} \cdot \nabla \right) \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} -\int_\Omega j \left(\hat{\mathbf{u}} \cdot \nabla \right) \overline{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x}  + \int_\Omega j \hat{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \hat{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega j \hat{\mathbf{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \hat{\mathbf{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s}
$$
__NOTE: In the linearized equations, the convective term is *not* integrated by parts.__

### Bilinear equations
The bilinear equations are required for sensitivity analyses \ldots

__NOTE: To be documented.__






## Compressible
Assumptions:
- primitive variables
- no additional source terms
- no gravitational forces

### Nonlinear equations
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
where $\mu_\textrm{eff}$ is the effective dynamic viscosity (see viscosity models for details) and where $\mathbf{I}$ is the identity tensor.

### Mean flow equations (to be reviewed by Simon)

We consider the flow field to be comprised of a  time-invariant base flow, which can be either a time-averaged flow or fixed point solution (base flow), and the perturbation, such that 

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\mathbf{q}'(\mathbf{x},t), \qquad \mathbf{q}=(\mathbf{u},p,\rho)^{T}
$$

Note that for the compressible equtions, the overbar represents a favre-abverage for the velocity and a Reynolds average for pressure and density.  

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





### Weak form
To solver the base flow equations with FELiCS we need to write it in weak form, reading
$$
\int_\Omega j \nabla \cdot \left(\mathbf{X}_\mathbf{u}^* \otimes \overline{\rho}\overline{\mathbf{u}} \right) \cdot \overline{\mathbf{u}} \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \overline{\rho}\left( \left(\overline{\mathbf{u}} \otimes \mathbf{X}_\mathbf{u}^*\right) \cdot \overline{\mathbf{u}} \right)\cdot \mathbf{n} \, \mathrm{d}\mathbf{s}  + \int_\Omega j \overline{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \overline{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega j \overline{\mathbf{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \overline{\mathbf{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s} = 0
$$
where the mean viscous stress tensor $\overline{\tau}$ is
$$
\mathbf{\overline{\tau}} = \overline{\mu}_\textrm{eff}[(\nabla + \nabla ^T)\overline{\mathbf{u}}]
$$
and where $\overline{\mu}_\textrm{eff}$ is the mean effective dynamic viscosity (see viscosity models for details).


### Linear equations
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
and where $\mu_\textrm{eff}'$ is the effective dynamic viscosity (see viscosity models for details).

#### Weak form
The weak form of the linearized Navier--Stokes equations with normal mode ansatz, as implemented in FELiCS, is
$$
\int_\Omega \omega \hat{\rho}X_p^* \mathrm{d}\mathbf{x} = \int_\Omega j \widehat{\rho\mathbf{u}} \cdot \nabla X_p^* \, \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \widehat{\rho\mathbf{u}} \cdot \mathbf{n} X_p^* \, \mathrm{d}\mathbf{s}
$$
and
$$
\int_\Omega \omega \overline{\rho} \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} = - \int_\Omega j \left(\overline{\rho}\overline{\mathbf{u}} \cdot \nabla \right) \hat{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} -\int_\Omega j \left(\overline{\rho}\hat{\mathbf{u}} \cdot \nabla \right) \overline{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} - \int_\Omega j \left(\hat{\rho}\overline{\mathbf{u}} \cdot \nabla \right) \overline{\mathbf{u}} \cdot \mathbf{X}_\mathbf{u}^* \, \mathrm{d}\mathbf{x} + \int_\Omega j \hat{p} \nabla \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} - \int_{\partial\Omega} j \hat{p} \mathbf{X}_\mathbf{u}^* \cdot \mathbf{n} \, \mathrm{d}\mathbf{s} - \int_\Omega j \hat{\mathbf{\tau}} : \nabla \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{x} + \int_{\partial\Omega} j \hat{\mathbf{\tau}} \mathbf{n} \cdot \mathbf{X}_\mathbf{u}^* \mathrm{d}\mathbf{s}
$$
with
$$
\widehat{\rho\mathbf{u}} = \overline{\rho}\mathbf{\hat{u}} + \hat{\rho}\mathbf{\overline{u}}
$$
__NOTE: In the linearized equations, the convective term is *not* integrated by parts.__
