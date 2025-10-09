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

where $\mathbf{q} = (\mathbf{u}, p, ρ, ...)^{T} \in \mathbb{R}^N $ is the state vector, $\mathcal{N}$ denotes the nonlinear operator of the system (e.g. Navier-Stokes), $\mathcal{B}$ is the weight operator.


We decompose the flow field into a time-invariant base flow, $\overline{\mathbf{q}} \in \mathbb{R}^N$, representing either a time-averaged state or a fixed-point solution, and a perturbation, $\mathbf{q}' \in \mathbb{R}^N$, such that

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\mathbf{q}'(\mathbf{x},t)
$$

Considering a stationary baseflow we arrive at a Linear Time-Invariant (LTI) dynamical system describing the perturbation:

$$
\mathcal{B}\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}')+\mathbf{f}',
$$

where

$$
\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\overline{\mathbf{q}}}\in \mathbb{R}^{N\times N}
$$
represents the nonlinear operator linearized around the about the base flow $\overline{\mathbf{q}}$ while the forcing $\mathbf{f}' = \mathcal{N}(\overline{\mathbf{q}})+O(|\mathbf{q'}|^2)+\mathbf{q} \in \mathbb{R}$ collects the nonlinear operator acting on the base flow and the nonlinear terms ([Rolandi et al. 2024](https://doi.org/10.1007/s00162-024-00717-x)). 



We consider harmonic forcing 
$$
\mathbf{f}' = \hat{\mathbf{f}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c., 
$$
and  harmonic responce
$$
\mathbf{q}' = \hat{\mathbf{q}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c. \ .
$$

Inserting in the LTI dynamical system and discretisation leads to the matrix formulation

$$
-\mathrm{j}\omega\mathbf{B} = \mathbf{A}\hat{\mathbf{q}}+\hat{\mathbf{f}}
$$

where $\mathbf{A}, \mathbf{B}$ are the discrete versions of the Jacobian $\mathcal{L}$ and limiter matrix $\mathcal{B}$, respectively.  

For resolvent analysis, this can be rearranged as

$$
\hat{\mathbf{q}} = \mathbf{R} \hat{\mathbf{f}},
$$

defining the resolvent operator

$$
\mathbf{R} = \left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1}.
$$

$\mathbf{R}$ acts as a transfer function that maps a given forcing $\hat{\mathbf{f}}$ to the linear response $\hat{\mathbf{q}}$.

## Optimal forcing-responce analysis

It is the goal to identify optimal input-output pairs $(\hat{\mathbf f}(\omega), \hat{\mathbf q}(\omega))$ maximizing amplification at frequency $\omega$.

This is formalized in the maximisation of the resolvent gains $\sigma$: 
$$
\sigma^2(\omega) = \max\limits_{\hat{\mathbf f}} \frac{\| \hat{\mathbf q}(\omega) \|^2}{\| \hat{\mathbf f}(\omega) \|^2}
$$ 

where $\|\cdot \|$ is a suitable energy norm detailed later.

We solve this optimisation problem by performing
the singular value decomposition (SVD) of the resolvent operator

$$\boxed{{\mathbf{R}(\omega) = Q \Sigma F^* } = \sum_j\hat{\mathbf q}_j\sigma_j\hat{\mathbf f}_j^*}.$$

This yields
|                           |                                       |                       |
|---------------------------|---------------------------------------|-----------------------|
| right singular vector:    | $F = \left[ \hat{\mathbf{f}}_1, \hat{\mathbf{f}}_2, ..., \hat{\mathbf{f}}_N \right]\ \in \mathbb{C}^{N\times N}$ | optimal forcing      |
| left singular vector:     | $Q = \left[ \hat{\mathbf{q}}_1, \hat{\mathbf{q}}_2, ..., \hat{\mathbf{q}}_N \right]\ \in \mathbb{C}^{N\times N}$ | optimal response     |
| singular values:          | $\Sigma = \mathrm{diag}(\sigma_1, \sigma_2, ..., \sigma_N )\ \in \mathbb{R}^{N\times N}$ | resolvent gain        |

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
\hat{q} = R W_{\text{FEM}} \hat{f}\ .
$$
We make the framework more flexible by introducing limiter operators ([see Towne et al. 2018](https://doi.org/10.1017/jfm.2018.283)), defining the projections:
$$ 
\hat{y} = P_r \hat{q}, \quad \hat{f} = P_f \hat{\eta},
$$
where $\hat{y}$ and $\hat{\eta}$ are the response and input of the reduced system, respectively. 
The operators $P_r$ and $P_f$ allow to select and weight the spatial regions and variables involved in the response and inputs, respectively.
The resolvent operator for the new system becomes: 
$$
\hat{y}=P_r \hat{q} = \tilde{R}\hat{\eta},\quad\text{with}\quad \mathbf{\tilde{R}=P_r R W_{\text{FEM}}P_f}$$




#### Inner product and energy norm

The inner product in the discretized domain is written as: 
$$
\langle \hat{a}, \hat{b} \rangle = \hat{a}^H W \hat{b},
$$
where $"\cdot^H"$ indicates the Hermitian transpose (transpose + complex conjugate).
The energy norm for the output term is then defined as $\|\hat{y}\|^2=\hat{y}^HW_{r}\hat{y}$ and the energy of the input term as $\|\hat{\eta}\|^2=\hat{\eta}^HW_{f}\hat{\eta}$. 


The gain squared is defined using this energy norm: 
$$
\sigma^2 = \frac{\|\hat{y}\|^2}{\| \hat{\eta}\|^2} = \frac{\hat{y}^H W_{r} \hat{y}}{\hat{\eta}^H W_{f} \hat{\eta}}.
$$

Currently, two norms are available by default in FELiCS:
- **TKE**: The turbulent kinetic energy norm is defined as $k = \int_{\Omega} \frac{1}{2}\bar{\rho}(u'^2+v'^2+w'^2)dx$, where $u',\,v',\,w'$ are the components of the fluctuation velocity vector. 
- **Chu**: The Chu norm for compressible fluctuations ([see George \& Sujith 2011](https://doi.org/10.1016/j.jsv.2011.06.016)), defined as $k=\int_{\Omega}\frac{1}{2}\bar{\rho}(u'^2+v'^2+w'^2)+\frac{1}{2}\frac{R\bar{T}}{\bar{\rho}}\rho'^2+\frac{1}{2}\frac{c_p\bar{\rho}}{\gamma\bar{T}}T'^2 dx$, where $R$ is the gas constant, $\gamma$ is the ratio of specific heats at constant pressure and volume, and $c_p$ is the specific heats at constant pressure.

The discretized form of these expressions are contained in the $W_f$ and $W_r$ operators, defining the norm for the forcing and response terms, respectively.

```{note}
The current version of the Chu norm implemented in FELiCS is only compatible with the state vector defined as $q=[\rho, u, T]$. 
```

### Equivalent EVP
We then introduce the resolvent in the expression for the gains $$\sigma^2=\frac{\hat{\eta}^H \tilde{R}^H W_{r} \tilde{R} \hat{\eta}}{\hat{\eta}^H W_{f} \hat{\eta}}.
$$
Using the Cholesky decomposition for the forcing weighting matrix,  $W_{f} =  M_f^H M_f$, and introducing a new function $\hat{g}=M_f\hat{\eta}$, we can re-write the definition of the gain as 
$$\sigma^2 = \frac{\hat{g}^H (M_f^H)^{-1} \tilde{R}^H W_{r} \tilde{R} (M_f)^{-1} \hat{g}}{\hat{g}^H\hat{g}}.$$
It takes on the form of a Rayleigh quotient: $\max \sigma^2$ is the solution of the following **HEVP (Hermitian Eigenvalue Problem)**:
$$
(M_f^H)^{-1} \tilde{R}^H W_{r} \tilde{R} (M_f)^{-1} \hat{g}=\lambda\hat{g},
$$
or in the original variables:
$$
W_{f}^{-1} \tilde{R}^H W_{r} \tilde{R} \hat{\eta}=\lambda\hat{\eta}$$
or:
$$W_{f}^{-1} P_f^H W_{\text{FEM}}^H R^H P_r^H W_{r}P_r R W_{\text{FEM}}P_f \hat{\eta}=\lambda\hat{\eta}.$$

```{note}
For real operators, such as $P_f,\, W_{FEM}, \,...$ the Hermitian transpose is just the transpose.
```

Finally, we can also re-write the H-EPV in terms of the linear operator:
$$
W_{f}^{-1}~P_f^H~ W_{\text{FEM}}^H~(\left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1})^H~P_r^H~W_{r}~P_r~\left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1}~W_{\text{FEM}}~P_f~\hat{\eta}=\lambda~\hat{\eta}.
$$
This is the expression implemented in FELiCS.

After solving the HEVP, the full forcing vectors are obtained from $\hat{f}=P_f\hat{\eta}$ and the response vectors are re-computed by applying $\hat{q} = R \hat{f}$.

```{note}
Based on the above, the forcing vectors $\hat{f}$ obtained from FELiCS have a unit norm over the (forcing) domain, while the response vectors $\hat{q}$ have a norm equal to $\sigma$.
```