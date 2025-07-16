# Input-output analysis

References:
- [Avdonin et al. 2018](https://doi.org/10.1016/j.proci.2018.06.142)
- [Kaiser et al. 2021](https://doi.org/10.1017/jfm.2021.151)
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Input-output analysis is a framework used to study how a linearized dynamical system responds to a given external forcing. Starting from the linearized equations around a steady base state, one models the system as an operator that maps input disturbances to system responses in frequency space. The focus is not on natural instabilities, but rather on how specific inputs (e.g. body forces, boundary perturbations) generate specific outputs (e.g. velocity or pressure fields). Unlike resolvent analysis, which typically uses singular value decomposition to identify optimal forcings, input-output analysis directly computes the response to a prescribed forcing using the resolvent operator. This approach is particularly useful when the structure or frequency content of the disturbance is known or imposed.

## Linearized Input-output System
We start with a general nonlinear equation, written in compact form as
$$
\frac{\mathrm{d}\mathbf{q}}{\mathrm{d}t}=\mathcal{N}(\mathbf{q})+\mathbf{C}\mathbf{f}
$$

- $\mathcal{N}$: nonlinear operator (e.g. navier-Stokes)
- $\mathbf{C}$: input matrix (maps external forcing $\mathbf{f}$ into the state space)
- $\mathbf{f}$: external forcing (e.g. actuation, volumne force)

Let $\overline{\mathbf{q}}$ be the steady base state , i.e. $\mathcal{N}(\overline{\mathbf{q}})=0$ and introduce a small pertubation 

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\epsilon\mathbf{q}'(\mathbf{x},t), \qquad \epsilon \ll 1.
$$

Linearizsation gives the linearized input-output system:

$$
\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}')+\mathbf{C}\mathbf{f}(t)
$$

where  $\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\overline{\mathbf{q}}}\in \mathbb{R}^{N\times N}$
is the Jacobian at the base state $\overline{\mathbf{q}}$. 

We further define the outpout equation 
$$\mathbf{y}' =\mathrm{D}\mathbf{q}'$$

with 
- $\mathbf{D}$: output matrix, projecting state to measurable quantities

### Harmonic forcing 
We assume the forcing and respoce to be periodic in time, reading

- harmonic forcing: $\mathbf{f}' = \hat{\mathbf{f}}\mathrm{e}^{-j\omega t} +c.c.$
- harmonic responce: $\mathbf{q}' = \hat{\mathbf{q}}\mathrm{e}^{-j\omega t} +c.c.$

Insertin in the linearized system of equations leads to

$$
-j\omega\mathbf{B}\hat{\mathbf{q}} = \mathbf{A}\hat{\mathbf{q}}+\mathbf{C}\hat{\mathbf{f}}.
$$


We solve for the output response  

$$
 \hat{\mathbf{q}} = (-j \omega \mathbf{B}-\mathbf{A})^{-1}\mathbf{C}\hat{\mathbf{f}}=\mathbf{R}(\omega)\mathbf{C}\hat{\mathbf{f}},
$$

with 

- resolvent operator:  $\mathbf{R} = (-j \omega \mathbf{B}-\mathbf{A})^{-1}$ .

The output is computed as  

$$
 \hat{\mathbf{y}}= \mathbf{D}\hat{\mathbf{q}} = \mathbf{D}(-j \omega \mathbf{B}-\mathbf{A})^{-1}\mathbf{C}\hat{\mathbf{f}}=\mathbf{H}(\omega)\mathbf{C}\hat{\mathbf{f}},
$$

where 

- input-output transfer functions:  $\mathbf{H} = \mathbf{D}\mathbf{R}\mathbf{C}$ .


## FELiCS implementation (to be done)




