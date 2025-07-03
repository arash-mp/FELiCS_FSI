# Implementation of steady incompressible RANS equations (with k-$\varepsilon$ turbulent model) - DONE (for future versions of FELiCS)

**Xiuyang has done a draft, and is reviewing. Welcome for any comments.**

## Nonlinear form of steady incompressible RANS equations
The starting point are the four equations (Continuity, Momentum, k and epsilon equations) for incompressible flow:
$$
    \nabla \cdot \mathbf{u} = 0
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
$$
where 
$$\nu_{Eff} = \nu +\nu_t$$
$$\nu_t = C_\mu \frac{k^2}{\varepsilon}$$
$$\nu_k = \nu + \frac{\nu_t}{\sigma_k}$$
$$P_k = \frac{1}{2}\nu_t[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]$$

$$\nu_\varepsilon = 
\nu+ \frac{\nu_t}{\sigma_\varepsilon}
$$
.

In the standard $k- \varepsilon$ model, the constants in the equations are listed as follow:
$C_\mu = 0.09$, $\sigma_k = 1.0$, $\sigma_\varepsilon = 1.3$, $C_{1\varepsilon } = 1.44$, $C_{2\varepsilon} = 1.92$.

## Linear form of steady incompressible RANS equations
To calculate base flow with Newton-Raphson method, or do linear analysis, it is necessary to get the linear form of the governing equations.
Therefore we derived the linear form of steady incompressible RANS equations and summerize them here.
The notation of $\delta \cdot$ (such as $\delta \mathbf{u}, \delta p, \delta k, \delta \varepsilon$, ...) presents the small perturbations on each variable.
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
where 
$$\delta \nu_{Eff} = \delta \nu_t = 2C_\mu \frac{k}{\varepsilon}\delta k - C_\mu \frac{k^2}{\varepsilon^2} \delta \varepsilon $$
$$\delta \nu_k = \frac{\delta \nu_t}{\sigma_k}$$
$$\delta P_k = \frac{1}{2}\delta \nu_t[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T] + \nu_t[\nabla\delta \mathbf{u} + {\nabla\delta\mathbf{u}}^T]:[\nabla\mathbf{u} + {\nabla\mathbf{u}}^T]$$
$$\delta \nu_\varepsilon = \frac{\delta \nu_t}{\sigma_\varepsilon}$$


## Variational form & Weak form Galerkin

To solve equation sets in FELiCS, it is necessary to write them in the following variational form
$$
\int (\mathcal{F} \cdot X) dx = 0
$$

Where $\mathcal{F}$ is the governing equation, $X$ is the Test function, their inner products are intergrated over the domain.

When calculating the integration expression (variational form), it is sometimes important to perform integration by part on some terms. 
This can not only decrease the order of spacial differences on specific parameters but also increase the rubostness of the solver. 
By doing so, the new expressions are usually called weak form Galerkin. 
This kind of equations are coded and implemented in FELiCS. 
The weak forms of linear and non-linear expressions for each equations are shown in the following.

### Weak form of continuity equation
The weak form of non-linear continuity equation is
$$
\int (\mathbf{u} \cdot \nabla X) dx - \int_\Omega (\mathbf{u} \cdot \mathbf{n}) X ds = 0
$$
where $\mathbf{n}$ is the unit vector which is normal to the boundary facets, $\int_\Omega (\cdot) ds$ is integration on boundaries.

The weak form for Linear continuity equation is
$$
\int(\delta \mathbf{u} \cdot \nabla X) dx - \int_\Omega (\delta \mathbf{u} \cdot \mathbf{n}) X ds = 0
$$

### Weak form of momentum equation
The weak form of NON-linear momentum equation is
$$
\int(
    \mathbf{u} \cdot (\nabla \cdot (\mathbf{X} \otimes \mathbf{u}) )
    + p \nabla \cdot (\mathbf{X})
    - \nu_{Eff}[\nabla\mathbf{u} + (\nabla\mathbf{u})^T] : \nabla\mathbf{X}
    + \frac{2}{3} k \nabla \cdot \mathbf{X}
) dx
\\
-\int_\Omega(
    \mathbf{u} \cdot (\mathbf{u} \otimes \mathbf{X}) \cdot \mathbf{n}
    + p \mathbf{X} \cdot \mathbf{n}
    - \nu_{Eff}[\nabla\mathbf{u} + (\nabla\mathbf{u})^T] : (\mathbf{X} \otimes\mathbf{n})
    + \frac{2}{3} k (\mathbf {X} \cdot \mathbf{n})
) ds
$$

The weak form for Linear momentum equation is 
$$
\int(
    \delta \mathbf{u} \cdot (\nabla \cdot (\mathbf{X} \otimes \mathbf{u}))
    + \mathbf{u} \cdot (\nabla \cdot (\mathbf{X} \otimes \delta \mathbf{u}))
    + \delta p (\nabla \cdot \mathbf{X})
    - \nu_{Eff} (\nabla \delta \mathbf{u} + (\nabla \delta \mathbf{u})^T) : \nabla \mathbf{X}
    - \delta \nu_{Eff} (\nabla \mathbf{u} + (\nabla \mathbf{u})^T) : \nabla \mathbf{X}
    + \frac{2}{3} \delta k (\nabla \cdot \mathbf{X})
)dx
\\
-\int_\Omega(
    ((\delta \mathbf{u} \otimes \mathbf{X}) \cdot \mathbf{u}) \cdot \mathbf{n}
    + ((\mathbf{u} \otimes \mathbf{X}) \cdot \delta \mathbf{u}) \cdot \mathbf{n}
    + \delta p \mathbf{X} \cdot \mathbf{n}
    - \nu_{Eff} (\nabla \delta \mathbf{u} + (\nabla \delta \mathbf{u})^T ) : (\mathbf{X} \otimes \mathbf{n})
    - \delta \nu_{Eff} (\nabla \mathbf{u} + (\nabla \mathbf{u})^T ) : (\mathbf{X} \otimes \mathbf{n})
    + \frac{2}{3} \delta k \mathbf{X} \cdot \mathbf{n}
)ds
$$

### Weak form of k equation

The weak form of NON-linear k equation is
$$
\int(
    (\nabla \cdot (\mathbf{u} X)) k
    - \nu_k (\nabla k \cdot \nabla X)
    + P_k X
    - \varepsilon X
)dx
\\
-\int_\Omega(
    (\mathbf{u} \cdot \mathbf{n}) X k
    - \nu_k (X (\nabla k \cdot \mathbf{n}))
)ds
$$

The weak form of Linear k equation is
$$
\int(
    (\nabla \cdot (\delta \mathbf{u}X))k
    +(\nabla \cdot (\mathbf{u} X)) \delta k
    -\nu_k (\nabla \delta k \cdot \nabla X)
    -\delta \nu_k (\nabla k \cdot \nabla X)
    + \delta P_k X
    -\delta \varepsilon X
)dx
\\
-\int_\Omega(
    ((\delta  \mathbf{u} X) \cdot \mathbf{n}) k 
    + ((\mathbf{u} X) \cdot \mathbf{n}) \delta k
    - \nu_k ((\nabla \delta k) \cdot \mathbf{n}) X
    - \delta \nu_k ((\nabla k) \cdot \mathbf{n}) X
)ds
$$

### Weak form of $\varepsilon$ equation

The weak form of non-linear $\varepsilon$ equation is

$$
\int(
    \varepsilon \nabla \cdot (\mathbf{u} X)
    - \nu_\varepsilon (\nabla \varepsilon \cdot \nabla X)
    + C_1 \frac{\varepsilon}{k} P_k X
    - C_2 \frac{\varepsilon ^2}{k} X
)dx
\\
-\int_\Omega(
    (\mathbf{u} \cdot \mathbf{n}) X \varepsilon
    -\nu_\varepsilon (\nabla \varepsilon \cdot \mathbf{n}) X
)ds
$$

The weak form of Linear $\varepsilon$ equation is
$$
\int(
    (\nabla \cdot (\delta \mathbf{u}X))\varepsilon
    +(\nabla \cdot (\mathbf{u} X)) \delta \varepsilon
    -\nu_\varepsilon (\nabla \delta \varepsilon \cdot \nabla X)
    -\delta \nu_\varepsilon (\nabla \varepsilon \cdot \nabla X)
    + C_1 \delta \varepsilon \frac{1}{k} P_k X
    - C_1 \delta k \frac{\varepsilon}{k^2} P_k X
    + C_1 \frac{\varepsilon}{k} \delta P_k X
    - 2 C_2 \delta \varepsilon \frac{\varepsilon}{k} X
    + C_2 \delta k \frac{\varepsilon ^2}{k^2} X
)dx
\\
-\int_\Omega(
    ((\delta  \mathbf{u} X) \cdot \mathbf{n}) \varepsilon 
    + ((\mathbf{u} X) \cdot \mathbf{n}) \delta \varepsilon
    - \nu_\varepsilon ((\nabla \delta \varepsilon) \cdot \mathbf{n}) X
    - \delta \nu_\varepsilon ((\nabla \varepsilon) \cdot \mathbf{n}) X
)ds
$$

