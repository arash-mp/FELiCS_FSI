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
- Mueller2024 JFM
- Barkley2006
- Towne2018
- Demange2024 AIAA


## Incompressible
Assumptions:
- $\rho = \textrm{const}$, this is reasonable if Mach number is low, no large temperature gradients/changes, no acoustics, no combustion
- no additional source terms
- no gravitational forces

### Nonlinear equations
$$
\nabla \cdot \mathbf{u} = 0
$$

### Linear equations
$$
\nabla \cdot \mathbf{u}' = 0
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

### Linear equations
$$
\frac{\partial \rho'}{\partial t} + \nabla \cdot (\overline{\rho} \mathbf{u}' + \rho' \overline{\mathbf{u}}) = 0
$$

### Bilinear equations