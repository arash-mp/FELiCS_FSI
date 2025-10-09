# Linear stability analysis
Linear stability analysis is a method used to determine whether small disturbances to a steady base state grow or decay over time. The system is linearized around the base flow, and solutions are sought in the form of exponentially growing or decaying modes. By solving an eigenvalue problem, one identifies the growth rates and shapes of these modes. If any mode grows over time, the base flow is considered unstable. This approach provides insight into the natural tendencies of the system to amplify disturbances without external forcing, focusing purely on the system’s internal dynamics.


**References:**
- [Barkley et al. 2006](https://doi.org/10.1209/epl/i2006-10168-7)
- [Sipp et al. 2010](https://doi.org/10.1115/1.4001478)
- [Kaiser et al. 2017](https://doi.org/10.1115/GT2017-63649)
- [Müller et al. 2020](https://doi.org/10.1017/jfm.2019.1063)



## Definition of the linear operator 
We start with a general nonlinear equation, written in compact form as
$$
\mathcal{B}\frac{\mathrm{d}\mathbf{q}}{\mathrm{d}t}=\mathcal{N}(\mathbf{q}),
$$

where $\mathbf{q} = (\mathbf{u}, p, ρ, ...)^{T} \in \mathbb{R}^N $ is the state vector (conservative variables), $\mathcal{N}$ denotes the nonlinear operator of the system (e.g. Navier-Stokes), $\mathcal{B}$ is a limiter operator; zero if the time derivative is not considered, else one.

We decompose the flow field into a time-invariant base flow, $\overline{\mathbf{q}} \in \mathbb{R}^N$, representing either a time-averaged state or a fixed-point solution, and a perturbation, $\mathbf{q}' \in \mathbb{R}^N$, such that 

$$
\mathbf{q}(\mathbf{x},t) = \overline{\mathbf{q}}(\mathbf{x})+\epsilon\mathbf{q}'(\mathbf{x},t)\ .
$$

In contrast to the resolvent analysis, we assume the  perturbation to be small, $\epsilon \ll 1$. 
Upon inserting this ansatz into the governing equations and linearisation, we arrive at the homogeneous equation for a Linear Time-Invariant (LTI) dynamical system describing the  perturbation:

$$
\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}'),
$$

where 

$$
\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\overline{\mathbf{q}}}\in \mathbb{R}^{N\times N}
$$
represents the Jacobian evaluated at the base state $\overline{\mathbf{q}}$. 


## Spectral analysis of the linear operator

We assume the perturbation takes the form of normal modes,
$$
\mathbf{q}'=\hat{\mathbf{q}}\mathrm{e^{-\mathrm{j}\omega t}}+c.c., 
$$
where $\hat{\mathbf{q}}$ denotes the complex amplitude, and the complex frequency is given by $\omega = \omega_r + \mathrm{j}\omega_i$.
Inserting this in the linearized perturbation equation leads to the eigenvalue problem,

$$
\mathbf{A} \hat{\mathbf{q}} = -\mathrm{j} \omega \mathbf{B} \hat{\mathbf{q}},
$$

where $\omega$ is the eigenvalue, $\hat{\mathbf{q}}$ is the eigenvector and $\mathbf{A}, \mathbf{B}$ are the discrete versions of the Jacobian $\mathcal{L}$ and $\mathcal{B}$, respectively.

Solving the eigenvalue problem yields complex eigenvalues, whose components have the following interpretations:

- $\omega_i > 0$: exponential temporal growth

- $\omega_i < 0$: exponential temporal decay

- $\omega_r$: oscillation frequency of the mode

- $\hat{\mathbf{q}}$: spatial structure (mode shape)


## FELiCS implementation 

The discretized eigenvalue problem is solved using a SLEPc-based solver implemented in FELiCS. The solver computes a set of eigenvalues $\omega_k$ ($k=1,2,...,n$) located nearest to a defined guess value $\omega_{\text{guess}}$. For each eigenvalue $\omega_k$, the corresponding eigenvector $\hat{\mathbf{q}}_k$ is obtained.

<!-- ### Derivation: For incompressible Navier-Stokes equations
Modal analysis can be used for linear stability analysis on flow governed by Navier-Stokes equations (NSE). By linearize NSE around a steady state (so called base flow), one can get the aformentioned linear system and solve it. Here we show the linearization on the incompressible NSE, which is:
$$\nabla \cdot \mathbf{u} = 0$$
$$\frac{d\mathbf{u}}{dt} + ((\mathbf{u} \cdot\nabla)\mathbf{u}) + \frac{1}{\rho} \nabla p - \nu \nabla^2 (\mathbf{u}) = 0$$
Perform double decompostion on the variables $\mathbf{u}$ and $p$,
$$\mathbf{u} = \bar{\mathbf{u}} + \mathbf{u}^\prime$$
$$p = \bar{p} + p^\prime$$
where $\bar{\mathbf{u}}$ and $\bar{p}$ are steady-state variables which can be get by solving the steady state equations (equaitons without time derivitive term). 
$\mathbf{u}^\prime$ and $p^\prime$ are time-dependent small fluctuations around the steady state. Considering the periodical characteristic of the fluctuations, they can be expressed like
$$\phi^\prime = \hat{\phi} e ^ {-j \omega t}$$
where $\hat{\phi}$ is the amplitude of fluctuation, and $\omega$ is the frequency.
Then the double decomposed $\mathbf{u}$ and $p$ can be expressed as:
$$\mathbf{u} = \bar{\mathbf{u}} + \hat{\mathbf{u}} e^{-j \omega t}$$
$$p = \bar{p} + \hat{p}e^{- j \omega t}$$
Then insert the decompostion into NSE and ignore the high-order fluctuations, the following linearized equations can be get:
$$(\nabla\cdot\hat{\mathbf{u}} )e^{-i \omega t} = 0$$
$$\hat{\mathbf{u}} $$

### Residual calculation
After calculation, the relative residuals for each pair of eigenvalues and eigenvectors are calculated as
$$
\mathcal{R}_i = \frac{||A\hat{q_i} - \omega_i B \hat{q_i}||}{||\omega_i B \hat{q_i}||} .
$$ -->
