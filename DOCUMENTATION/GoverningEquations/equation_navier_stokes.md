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


## Incompressible
Assumptions:
- $\rho = \textrm{const}$, this is reasonable if Mach number is low, no large temperature gradients/changes, no acoustics, no combustion
- no additional source terms
- no gravitational forces

### Nonlinear equations
$$
\nabla \cdot \mathbf{u} = 0
$$

$$
\rho \frac{\mathbf{u}}{\partial t} + \rho(\mathbf{u}\cdot \nabla) \mathbf{u} + \nabla p - \nabla \cdot \tau = 0
$$

where the viscous stress tensor $\tau$ is
$$
\tau = \mu(\nabla + \nabla ^T)\mathbf{u}
$$
### Linear equations
$$
\nabla \cdot \mathbf{u}' = 0
$$

$$
\rho \frac{\partial \mathbf{u}'}{\partial t} + \rho [(\overline{\mathbf{u}}\cdot \nabla) \mathbf{u}' + (\mathbf{u}' \cdot \nabla) \overline{\mathbf{u}}] + \nabla p' - \nabla \cdot \tau ' = 0
$$

where the fluctuated viscous stress tensor $\tau '$ is
$$
\tau' = \mu'[(\nabla + \nabla ^T)\overline{\mathbf{u}}]+ \overline{\mu}[(\nabla + \nabla ^T) \mathbf{u}']
$$

### Bilinear equations



## Compressible
Assumptions:
- no additional source terms
- no gravitational forces

### Nonlinear equations
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
$$

$$
\frac{\rho\mathbf{u}}{\partial t} + \rho(\mathbf{u}\cdot \nabla) \mathbf{u} + \nabla p - \nabla \cdot \tau = 0
$$

where the viscous stress tensor $\tau$ is
$$
\tau = \mu[(\nabla + \nabla ^T)\mathbf{u} - \frac{2}{3} (\nabla \cdot \mathbf{u})\mathbf{I}]
$$
where $\mathbf{I}$ is the identity tensor.


### Linear equations
$$
\frac{\partial \rho'}{\partial t} + \nabla \cdot (\overline{\rho} \mathbf{u}' + \rho' \overline{\mathbf{u}}) = 0
$$

$$
\overline{\rho} \frac{\partial \mathbf{u}'}{\partial t} + \overline{\rho} [(\overline{\mathbf{u}}\cdot \nabla) \mathbf{u}' + (\mathbf{u}' \cdot \nabla) \overline{\mathbf{u}}] + \rho' (\overline{\mathbf{u}}\cdot \nabla)\overline{\mathbf{u}} + \nabla p' - \nabla \cdot \tau ' = 0
$$
where the fluctuated viscous stress tensor $\tau'$ is

$$
\tau' = \mu'[(\nabla + \nabla ^T)\overline{\mathbf{u}} - \frac{2}{3} \nabla \cdot \overline{\mathbf{u}} \mathbf{I}]+ \overline{\mu}[(\nabla + \nabla ^T) \mathbf{u}' - \frac{2}{3}( \nabla \cdot \mathbf{u}' )\mathbf{I}]
$$



### Bilinear equations
The bilinear forms for COMPRESSIBLE Navier--stokes equations are not implemented yet.