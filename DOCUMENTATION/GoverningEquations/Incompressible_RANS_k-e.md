# Implementation of steady incompressible RANS equations (with k-$\varepsilon$ turbulent model)

## Xiuyang is working on this documentation, and it is not done yet.

## Nonlinear form of steady incompressible RANS equations
The starting point are the four equations (Continuity, Momentum, k and epsilon equations) for incompressible flow:
$$
    \nabla \mathbf{u} = 0
$$
$$
    (\mathbf{u} \cdot \nabla) \mathbf{u}
    =
    - \frac{1}{\rho} \nabla p
    + \nabla \cdot \nu_{Eff} [\nabla \mathbf{u} + (\nabla \mathbf{u})^T]
    - \frac{2}{3} \nabla k
$$
$$
    \mathbf{u} \cdot \nabla k
    =
    \nabla \cdot (\nu_k \nabla k)
    + P_k
    - \varepsilon
$$
$$
    \mathbf{u} \cdot \nabla \varepsilon
    =
    \nabla \cdot (\nu_\varepsilon \nabla \varepsilon)
    + C_{1\varepsilon} \frac{\varepsilon}{k} P_k
    - C_{2\varepsilon} \frac{\varepsilon^2}{k}
    - S_{\varepsilon}
$$
where $\nu_{Eff} = \nu +\nu_t$, 
$\nu_t = C_\mu \frac{k^2}{\varepsilon}$,
$\nu_k = \nu + \frac{\nu_t}{\sigma_k}$, $P_k = \frac{1}{2}\nu_t[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]$, 
$\nu_\varepsilon = 
\nu+ \frac{\nu_t}{\sigma_\varepsilon}
$.

In the standard $k- \varepsilon$ model, the constants in the equations are listed as follow:
$C_\mu = 0.09$, $\sigma_k = 1.0$, $\sigma_\varepsilon = 1.3$, $C_{1\varepsilon } = 1.44$, $C_{2\varepsilon} = 1.92$.

## Linear form of steady incompressible RANS equations
To calculate base flow with Newton-Raphson method, or do linear analysis, it is necessary to get the linear form of the governing equations.
Therefore we derived the linear form of steady incompressible RANS equations and summerize them here.
The notation of $\delta \cdot$ ($\delta \mathbf{u}, \delta p, \delta k, \delta \varepsilon$, ...) presents the small perturbations of each variable.
$$
    \nabla \delta \mathbf{u} = 0
$$
$$
    (\delta \mathbf{u} \cdot \nabla) \mathbf{u}
    + (\mathbf{u} \cdot \nabla) \delta \mathbf{u}
    =
    - \frac{1}{\rho} \nabla \delta p
    + \nabla \cdot \delta \nu_{Eff} [\nabla \mathbf{u} + (\nabla \mathbf{u})^T]
    + \nabla \cdot \nu_{Eff} [\nabla \delta \mathbf{u} + (\nabla \delta \mathbf{u})^T]
    - \frac{2}{3} \nabla \delta k
$$
$$
    \delta \mathbf{u} \cdot \nabla k
    + \mathbf{u} \cdot \nabla \delta k
    =
    \nabla \cdot (\delta \nu_k \nabla k)
    + \nabla \cdot (\nu_k \nabla \delta k)
    + \delta P_k
    - \delta \varepsilon
$$
$$
    \delta \mathbf{u} \cdot \nabla \varepsilon
    + \mathbf{u} \cdot \nabla \delta \varepsilon
    =
    \nabla \cdot (\delta \nu_\varepsilon \nabla \varepsilon)
    + \nabla \cdot (\nu_\varepsilon \nabla \delta \varepsilon)
    + C_{1\varepsilon} \frac{\delta \varepsilon}{k} P_k
    - C_{1\varepsilon} \frac{\varepsilon}{k^2} \delta k P_k
    + C_{1\varepsilon} \frac{\varepsilon}{k} \delta P_k
    - 2 C_{2\varepsilon} \frac{\varepsilon}{k} \delta \varepsilon
    + C_{2\varepsilon} \frac{\varepsilon^2}{k^2} \delta k
$$
where $\delta \nu_{Eff} = \delta \nu_t = 2C_\mu \frac{k}{\varepsilon}\delta k - C_\mu \frac{k^2}{\varepsilon^2} \delta \varepsilon $, 

$\delta \nu_k = \frac{\delta \nu_t}{\sigma_k}$, 

$\delta P_k = \frac{1}{2}\delta \nu_t[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T] + \nu_t[\nabla\delta \mathbf{u} + {\nabla\delta\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]$, 

$\delta \nu_\varepsilon = \frac{\delta \nu_t}{\sigma_\varepsilon}$.


!!! IMPORTANT TO REMEMBER: When write the weak form, remember to keep consisttency with FELiCS codes. Especially the convecting terms.

