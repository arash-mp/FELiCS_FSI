# Implementation of steady incompressible RANS equations (with k-$\varepsilon$ turbulent model)

## Nonlinear form of steady incompressible RANS equations
The starting point are the four equations (Continuity, Momentum, k and epsilon equations) in incompressible form:
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
$, the constants in the equations are listed as follow:

!!! IMPORTANT TO REMEMBER: When write the weak form, remember to keep consisttency with FELiCS codes. Especially the convecting terms.


## Resolvent Analysis for a Linear System 
The starting point is the equation for a linear system: 
$$
(A + j \omega B) \hat{q} = J\hat{q} = \hat{f},
$$
where $A$ is the linear operator. Applying a discretization scheme and considering a finite element method (FEM) weighting:  
$$
\hat{q} = R W_{\text{FEM}} \hat{f},\quad \text{with}\quad \mathbf{R=J^{-1}}.
$$
We make the framework more flexible by introducing limiter operators (similar to @towneSpectralProperOrthogonal2018), defining the projections:
$$ 
\hat{y} = P_r \hat{q}, \quad \hat{f} = P_f \hat{\eta},
$$
where $\hat{y}$ and $\hat{\eta}$ are the response and input of the reduced system, respectively. 
The $P_r$ and $P_f$ allow to select the spatial regions and variables involved in the response and inputs, respectively.
The resolvent operator for the new system then becomes: $$\hat{y}=P_r \hat{q} = \tilde{R}\hat{\eta},\quad\text{with}\quad \mathbf{\tilde{R}=P_r R W_{\text{FEM}}P_f}$$
## Inner Product and Energy Norm
The inner product in the discretized domain is written as: $$\langle \hat{a}, \hat{b} \rangle = \hat{a}^H W \hat{b},$$where $"\cdot^H"$ indicates the Hermitian transpose (transpose + complex conjugate).
The energy norm for the output term is then defined as $\|\hat{y}\|=\hat{y}^HW_{r}\hat{y}$ for the output term. 

## Definition of the Gain 
The gain squared is defined as: $$\sigma^2 = \frac{\|\hat{y}\|^2}{\| \hat{\eta}\|^2} = \frac{\hat{y}^H W_{r} \hat{y}}{\hat{\eta}^H W_{f} \hat{\eta}}$$
## SDV and equaivalent EVP
The we introduce the resolvent in the expression for the gains $$\sigma^2=\frac{\hat{\eta}^H \tilde{R}^H W_{r} R_2 \hat{\eta}}{\hat{\eta}^H W_{f} \hat{\eta}}$$
Using the Cholesky decomposition for the forcing weighting matrix,  $W_{in} =  M_f^H M_f$, and introducing a new function $\hat{g}=M_f\hat{\eta}$, we can re-write the definition of the gain as $$\sigma^2 = \frac{\hat{g}^H (M_f^H)^{-1} \tilde{R}^H W_{r} \tilde{R} (M_f)^{-1} \hat{g}}{\hat{g}^H\hat{g}},$$which takes on the form of a Rayleigh quotient: $\max \sigma^2$ is then the solution of the following **HEVP (Hermitian Eigenvalue Problem)**: $$(M_f^H)^{-1} \tilde{R}^H W_{r} \tilde{R} (M_f)^{-1} \hat{g}=\lambda\hat{g},$$or in the original variables:$$(W_{f})^{-1} \tilde{R}^H W_{r} \tilde{R} \hat{\eta}=\lambda\hat{\eta}$$or:$$(W_{f})^{-1} P_f^H W_{\text{FEM}}^H R^H P_r^H W_{r} P_r R W_{\text{FEM}}P_f \hat{\eta}=\lambda\hat{\eta}$$
## Final Expression

Finally, we can also re-write the H-EPV in terms of the linear operator, which is the expression we implement in FELiCS: $$(W_{f})^{-1}~P_f^H~ W_{\text{FEM}}^H~(J^{-1})^H~P_r^H~W_{r}~P_r~J^{-1}~W_{\text{FEM}}~P_f~\hat{\eta}=\lambda~\hat{\eta}$$
We can then get the full forcing from $\hat{f}=P_f\hat{\eta}$ and multiply it by the original resolvent operator $R$ to recover the full quantity vector $\hat{q}$. 
Note: for real operators, such as $P_f,\,W_{FEM},\,...$ the Hermitian transpose is just the transpose.

## Dimensions of the operators

Setting $N$ the number of degrees of freedom of the linear operator:
* The state and forcing vectors, $\hat{q}$ and $\hat{f}$, are of length $N$.
* The linear operator and thus the initial resolvent operator $R$ are square matrices of dimension $[N_1\times N_1]$.
* Assuming the input $\hat{\eta}$ and output $\hat{y}$ to be of length $N_{\eta}$ and $N_{y}$ respectively (they do not necessarily have the same dimension).
* Then, the dimensions of the limitor operators must be $P_r:[N_y \times N_1]$ and $P_f:[N_1 \times N_{\eta}]$. **NOTE: This is weird because in the old FELiCS the limitor operator for the response is a square matrix.** Might be combined with another operator for it to be the case (The forcing limitor operator has the right size in FELiCS).
* The weighting operators used to define the input and output norms must necessary be square matrices of dimensions $W_{in}:[N_{\eta} \times N_{\eta}]$ and $W_{out}:[N_y \times N_{y}]$

