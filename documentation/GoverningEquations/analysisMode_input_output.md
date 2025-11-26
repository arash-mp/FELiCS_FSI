# Input-output analysis

References:
- [Avdonin et al. 2018](https://doi.org/10.1016/j.proci.2018.06.142)
- [Kaiser et al. 2021](https://doi.org/10.1017/jfm.2021.151)
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Input-output analysis is a framework used to study how a linearized dynamical system responds to a given external forcing. Starting from the linearized equations around a steady base state, we model the system as an operator that maps input disturbances to system responses in frequency domain. The focus is not on natural instabilities, but rather on how specific inputs (e.g. body forces, boundary perturbations) generate specific outputs (e.g. velocity or pressure fields). Unlike resolvent analysis, which typically uses singular value decomposition to identify optimal forcings, input-output analysis directly computes the response to a prescribed forcing using the resolvent operator. This approach is particularly useful when the structure or frequency content of the disturbance is known or imposed.

## Linearized input-output system
The input-output analysis is analogous to the [resolvent analysis](./analysisMode_resolvent.md), except that the periodic forcing, $\mathbf{f}'$, is known a priori and chosen by the user. The linearized system that links the periodic forcing to the response of the flow in terms of fluctuations $\mathbf{q}'$ can be expressed as

$$
\mathcal{B}\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}')+\mathbf{f}',
$$

with the restriction matrix $\mathcal{B}$ and the linear operator $\mathcal{L}$ defined as

$$
\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\overline{\mathbf{q}}}\in \mathbb{R}^{N\times N}.
$$

$\mathcal{L}$ represents the Jacobian of the nonlinear operator linearized about the base flow $\overline{\mathbf{q}}$ while the forcing, unlike in the resolvent analysis, can be either a volume forcing or a forcing at a domain boundary. 

We further consider a harmonic forcing 
$$
\mathbf{f}' = \hat{\mathbf{f}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c., 
$$
and a harmonic response
$$
\mathbf{q}' = \hat{\mathbf{q}}\mathrm{e}^{-\mathrm{j}\omega t} +c.c. \ .
$$

Inserting in the linearized dynamical system and discretization leads to the matrix formulation

$$
-\mathrm{j}\omega\mathbf{B} = \mathbf{A}\hat{\mathbf{q}}+\hat{\mathbf{f}}
$$

where $\mathbf{A}, \mathbf{B}$ are the discrete versions of the Jacobian $\mathcal{L}$ and restriction matrix $\mathcal{B}$, respectively.  

Analogous to the resolvent analysis, this can be rearranged as

$$
\hat{\mathbf{q}} = \mathbf{R} \hat{\mathbf{f}},
$$

with the resolvent operator

$$
\mathbf{R} = \left(-\mathrm{j} \omega \mathbf{B}-\mathbf{A}\right)^{-1}.
$$

$\mathbf{R}$ acts as a transfer function that maps a given forcing $\hat{\mathbf{f}}$ to the linear response $\hat{\mathbf{q}}$.

## FELiCS implementation

The resolvent operator maps a mass-weighted forcing field $\hat{\mathbf{f}}^*$ to the corresponding linear response  $\hat{\mathbf{q}}$ of the flow via

$$
\hat{\mathbf{q}} = \mathbf{R} \hat{\mathbf{f}}^*\ .
$$

The mass weighted forcing, $\hat{\mathbf{f}}^*$, represents either a volume force or 
a boundary forcing. For a volume (body) force, the user-specified forcing $\hat{\mathbf{f}}$ is 
projected onto the finite-element space by multiplication with the mass matrix,

$$
\hat{\mathbf{f}}^* =\mathbf{W}_{\text{FEM}} \hat{\mathbf{f}}\ .
$$

For a boundary forcing, the weighted forcing is obtained directly through the weak 
imposition of the forced boundary values using Nitsche’s method. In practice, this is 
equivalent to supplying the specified forcing as a boundary condition to the 
linearized equations.

```{note}
To allow for the weak imposition, the boundary conditions for the forced quantities need to be set to "None" at the respective boundaries.
```

In both cases—volume forcing or boundary forcing—the forcing field is interpreted 
through the mean (base) flow, as discussed for example in [Tutorial 4](../Tutorials/input_output_analysis.md). 





