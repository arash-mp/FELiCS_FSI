# Resolvent analysis

Resolvent analysis is a linear systems approach used to understand how a dynamical system responds to harmonic forcing. The method starts by linearizing the governing equations around a steady base state, then analyzing how the system amplifies input disturbances at different frequencies. By treating the linearized operator as a transfer function, one can identify which inputs (forcings) lead to the strongest outputs (responses). This is done using a singular value decomposition of the resolvent operator, which reveals the most amplified structures and quantifies the gain. The method captures both modal and nonmodal amplification, making it especially powerful for studying flows with strong non-normal behavior.

**References:**
- [Beneddine et al. 2016](https://doi.org/10.1017/jfm.2016.331)
- [Towne et al. 2018](https://doi.org/10.1017/jfm.2018.675)
- [von Saldern et al. 2024](https://doi.org/10.1017/jfm.2024.922)
- [Müller et al. 2024](https://doi.org/10.1017/jfm.2024.679)



## Definition of the Resolvent Operator

We start with a general nonlinear equation, written in compact form as
$$
\mathcal{B}\frac{\mathrm{d}\mathbf{q}}{\mathrm{d}t}=\mathcal{N}(\mathbf{q}),
$$

where $\mathbf{q} = (\mathbf{u}, p, ρ, ...)^{T} \in \mathbb{R}^N $ is the state vector, $\mathcal{N}$ denotes the nonlinear operator of the system (e.g. Navier-Stokes), $\mathcal{B}$ is a restriction operator.


We decompose the flow field into a time-invariant base flow, $\overline{\mathbf{q}} \in \mathbb{R}^N$, representing either a time-averaged state or a fixed-point solution, and a perturbation, $\mathbf{q}' \in \mathbb{R}^N$, such that

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\mathbf{q}'(\mathbf{x},t)
$$

Considering a stationary base flow we arrive at a linear time-invariant dynamical system describing the perturbation:

$$
\mathcal{B}\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}')+\mathbf{f}',
$$

where

$$
\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\overline{\mathbf{q}}}\in \mathbb{R}^{N\times N}
$$
represents the nonlinear operator linearized around the about the base flow $\overline{\mathbf{q}}$ while the forcing $\mathbf{f}' = \mathcal{N}(\overline{\mathbf{q}})+O(|\mathbf{q'}|^2)+\mathbf{g} \in \mathbb{R}^N$ collects the nonlinear operator acting on the base flow and the nonlinear terms as well as an external forcing $\mathbf{g}$ ([Rolandi et al. 2024](https://doi.org/10.1007/s00162-024-00717-x)). 



We consider a harmonic forcing 
$$
\mathbf{f}' = \hat{\mathbf{f}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c., 
$$
and a harmonic responce
$$
\mathbf{q}' = \hat{\mathbf{q}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c. \ .
$$

Inserting in the linear time-invariant dynamical system and discretization leads to the matrix formulation

$$
-\mathrm{j}\omega\mathbf{B} = \mathbf{A}\hat{\mathbf{q}}+\hat{\mathbf{f}}
$$

where $\mathbf{A}, \mathbf{B}$ are the discrete versions of the Jacobian $\mathcal{L}$ and restriction matrix $\mathcal{B}$, respectively.  

For resolvent analysis, this can be rearranged as

$$
\hat{\mathbf{q}} = \mathbf{R} \hat{\mathbf{f}},
$$

defining the resolvent operator

$$
\mathbf{R} = \left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1}.
$$

$\mathbf{R}$ acts as a transfer function that maps a given forcing $\hat{\mathbf{f}}$ to the linear response $\hat{\mathbf{q}}$.

## Optimal forcing-response analysis

In the resolvent analysis, also optimal forcing-response analysis, the goal is to identify optimal input-output pairs $(\hat{\mathbf f}(\omega), \hat{\mathbf q}(\omega))$ maximizing amplification at a given frequency $\omega$.

This is formalized in the maximization of the resolvent gains $\sigma$: 
$$
\sigma^2(\omega) = \max\limits_{\hat{\mathbf f}} \frac{\| \hat{\mathbf q}(\omega) \|^2}{\| \hat{\mathbf f}(\omega) \|^2}
$$ 

where $\|\cdot \|$ is a suitable energy norm detailed later.

We solve this optimisation problem by performing
the singular value decomposition (SVD) of the resolvent operator

$$\boxed{{\mathbf{R}(\omega) = \mathbf{Q} \mathbf{\Sigma} \mathbf{F}^* } = \sum_j\hat{\mathbf q}_j\sigma_j\hat{\mathbf f}_j^*}.$$

This yields
|                           |                                       |                       |
|---------------------------|---------------------------------------|-----------------------|
| right singular vector:    | $\mathbf{F} = \left[ \hat{\mathbf{f}}_1, \hat{\mathbf{f}}_2, ..., \hat{\mathbf{f}}_N \right]\ \in \mathbb{C}^{N\times N}$ | optimal forcing      |
| left singular vector:     | $\mathbf{Q} = \left[ \hat{\mathbf{q}}_1, \hat{\mathbf{q}}_2, ..., \hat{\mathbf{q}}_N \right]\ \in \mathbb{C}^{N\times N}$ | optimal response     |
| singular values:          | $\mathbf{\Sigma} = \mathrm{diag}(\sigma_1, \sigma_2, ..., \sigma_N )\ \in \mathbb{R}^{N\times N}$ | resolvent gain        |

```{note}
The gain is sorted by decreasing order $\sigma_1\geq\sigma_2\geq ... \geq \sigma_N\geq 0$, with the resolvent spectral norm $\| \textbf{R} \| =\sigma_1$. 
```

```{note}
In FELiCS, the values returned in the 'gains.csv' file are the gain squared: $\sigma^2$. 
```


## FELiCS implementation 
In FELiCS, the resolvent implementation contains some additional utilities. They stem from more detailed definitions of the response and forcing norms. For example, it is possible to define spatial regions, in which the norm should be computed ("spatial restrictors") or weight spatial regions differently. In the following, it is described how those limiters are oncorporated in the resolvent formulation.


#### Additional weighting and limiter operators

Applying a discretization scheme and considering a finite element method weighting:  
$$
\hat{\mathbf{q}} = \mathbf{R} \mathbf{W}_{\text{FEM}} \hat{\mathbf{f}}\ .
$$
We make the framework more flexible by introducing limiter operators ([see Towne et al. 2018](https://doi.org/10.1017/jfm.2018.283)), defining the projections:
$$ 
\hat{\mathbf{y}} = \mathbf{P}_r \hat{\mathbf{q}}, \quad \hat{\mathbf{f}} = \mathbf{P}_f \hat{\boldsymbol{\eta}},
$$
where $\hat{\mathbf{y}}$ and $\hat{\boldsymbol{\eta}}$ are the response and input of the reduced system, respectively. 
The operators $\mathbf{P}_r$ and $\mathbf{P}_f$ allow to select and weight the spatial regions and variables involved in the response and inputs, respectively.
The resolvent operator for the new system becomes: 
$$
\hat{\mathbf{y}}=\mathbf{P}_r \hat{\mathbf{q}} = \tilde{\mathbf{R}}\hat{\boldsymbol{\eta}},\quad\text{with}\quad \mathbf{\tilde{\mathbf{R}}=\mathbf{P}_r R \mathbf{W}_{\text{FEM}}\mathbf{P}_f}$$




#### Inner product and energy norm

The inner product in the discretized domain is written as: 
$$
\langle \hat{\mathbf{a}}, \hat{\mathbf{b}} \rangle = \hat{\mathbf{a}}^H \mathbf{W} \hat{\mathbf{b}},
$$
where $"\cdot^H"$ indicates the Hermitian transpose (transpose + complex conjugate).
The energy norm for the output term is then defined as $\|\hat{\mathbf{y}}\|^2=\hat{\mathbf{y}}^H\mathbf{W}_r\hat{\mathbf{y}}$ and the energy of the input term as $\|\hat{\boldsymbol{\eta}}\|^2=\hat{\boldsymbol{\eta}}^H\mathbf{W}_f\hat{\boldsymbol{\eta}}$. 


The gain squared is defined using this energy norm: 
$$
\sigma^2 = \frac{\|\hat{\mathbf{y}}\|^2}{\| \hat{\boldsymbol{\eta}}\|^2} = \frac{\hat{\mathbf{y}}^H \mathbf{W}_r \hat{\mathbf{y}}}{\hat{\boldsymbol{\eta}}^H \mathbf{W}_f \hat{\boldsymbol{\eta}}}.
$$

Currently, two norms are available by default in FELiCS:
- **TKE**: The turbulent kinetic energy norm is defined as $k = \int_{\Omega} \frac{1}{2}\bar{\rho}(u'^2+v'^2+w'^2)dx$, where $u',\,v',\,w'$ are the components of the fluctuation velocity vector. 
- **Chu**: The Chu norm for compressible fluctuations ([see George \& Sujith 2011](https://doi.org/10.1016/j.jsv.2011.06.016)), defined as $k=\int_{\Omega}\frac{1}{2}\bar{\rho}(u'^2+v'^2+w'^2)+\frac{1}{2}\frac{R\bar{T}}{\bar{\rho}}\rho'^2+\frac{1}{2}\frac{c_p\bar{\rho}}{\gamma\bar{T}}T'^2 dx$, where $R$ is the gas constant, $\gamma$ is the ratio of specific heats at constant pressure and volume, and $c_p$ is the specific heats at constant pressure.

The discretized form of these expressions are contained in the $\mathbf{W}_f$ and $\mathbf{W}_r$ operators, defining the norm for the forcing and response terms, respectively.

```{note}
The current version of the Chu norm implemented in FELiCS is only compatible with the state vector defined as $q=[\rho, u, T]$. 
```

### Equivalent eigenvalue problem
We introduce the resolvent in the expression for the gains
$$
\sigma^2=\frac{\hat{\boldsymbol{\eta}}^H \tilde{\mathbf{R}}^H \mathbf{W}_r \tilde{\mathbf{R}} \hat{\boldsymbol{\eta}}}{\hat{\boldsymbol{\eta}}^H \mathbf{W}_f \hat{\boldsymbol{\eta}}}.
$$
Using the Cholesky decomposition for the forcing weighting matrix,  $\mathbf{W}_f =  \mathbf{M}_f^H \mathbf{M}_f$, and introducing a new function $\hat{\mathbf{g}}=\mathbf{M}_f\hat{\boldsymbol{\eta}}$, we can re-write the definition of the gain as 
$$\sigma^2 = \frac{\hat{\mathbf{g}}^H (\mathbf{M}_f^H)^{-1} \tilde{\mathbf{R}}^H \mathbf{W}_r \tilde{\mathbf{R}} (\mathbf{M}_f)^{-1} \hat{\mathbf{g}}}{\hat{\mathbf{g}}^H\hat{\mathbf{g}}}.$$
It takes on the form of a Rayleigh quotient: $\max \sigma^2$ is the solution of the following HEVP (Hermitian Eigenvalue Problem):
$$
(\mathbf{M}_f^H)^{-1} \tilde{\mathbf{R}}^H \mathbf{W}_r \tilde{\mathbf{R}} (\mathbf{M}_f)^{-1} \hat{\mathbf{g}}=\lambda\hat{\mathbf{g}},
$$
or in the original variables:
$$
\mathbf{W}_f^{-1} \tilde{\mathbf{R}}^H \mathbf{W}_r \tilde{\mathbf{R}} \hat{\boldsymbol{\eta}}=\lambda\hat{\boldsymbol{\eta}}$$
or:
$$\mathbf{W}_f^{-1} \mathbf{P}_f^H \mathbf{W}_{\text{FEM}}^H \mathbf{R}^H \mathbf{P}_r^H \mathbf{W}_r\mathbf{P}_r \mathbf{R} \mathbf{W}_{\text{FEM}}\mathbf{P}_f \hat{\boldsymbol{\eta}}=\lambda\hat{\boldsymbol{\eta}}.$$

```{note}
For real operators, such as $\mathbf{P}_f,\, W_{FEM}, \,...$ the Hermitian transpose is just the transpose.
```

Finally, we can also re-write the HEVP in terms of the linear operator:
$$
\mathbf{W}_f^{-1}~\mathbf{P}_f^H~ \mathbf{W}_{\text{FEM}}^H~(\left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1})^H~\mathbf{P}_r^H~\mathbf{W}_r~\mathbf{P}_r~\left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1}~\mathbf{W}_{\text{FEM}}~\mathbf{P}_f~\hat{\boldsymbol{\eta}}=\lambda~\hat{\boldsymbol{\eta}}.
$$
This is the expression implemented in FELiCS.

After solving the HEVP, the full forcing vectors are obtained from $\hat{\mathbf{f}}=\mathbf{P}_f\hat{\boldsymbol{\eta}}$ and the response vectors are re-computed by applying $\hat{\mathbf{q}} = \mathbf{R} \hat{\mathbf{f}}$.

```{note}
Based on the above, the forcing vectors $\hat{\mathbf{f}}$ obtained from FELiCS have a unit norm over the (forcing) domain, while the response vectors $\hat{\mathbf{q}}$ have a norm equal to $\sigma$.
```