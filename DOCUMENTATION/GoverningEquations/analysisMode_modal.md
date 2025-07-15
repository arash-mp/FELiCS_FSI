# Linear Stability Analysis


Linear stability analysis is a method used to determine whether small disturbances to a steady base state grow or decay over time. The system is linearized around the base flow, and solutions are sought in the form of exponentially growing or decaying modes. By solving an eigenvalue problem, one identifies the growth rates and shapes of these modes. If any mode grows over time, the base flow is considered unstable. This approach provides insight into the natural tendencies of the system to amplify disturbances without external forcing, focusing purely on the system’s internal dynamics.


## Linear Operator 
We start with a general nonlinear equation, written in compact form as
$$
\frac{\mathrm{d}\mathbf{q}}{\mathrm{d}t}=\mathcal{N}(\mathbf{q})+\mathbf{g}
$$

where $ \mathbf{q} = (\mathbf{u}, p, ρ, ...)^{T} \in  \mathbb{R}^N $ represents the state vector (conservative variables) and $\mathcal{N} demontes the nonlinear operator (e.g. the NS evolution operator).

We consider the flow field to be comprised of a  time-invariant base flow $\mathbf{q}_b \in \mathbb{R}^N  $, which can be either a time-averaged flow or fixed point solution, and the perturbation $ \mathbf{q}' \in \mathbb{R} $, such that 

$$
\mathbf{q}(\mathbf{x},t) = \mathbf{q}_b(\mathbf{x})+\epsilon\mathbf{q}'(\mathbf{x},t)
$$

In contrast to the resolvent analysis, we assume the pertubation to me small, $\epsilon \ll 1$. 
Upon inserting this ansatz into the governing equations and linearisation, we  arrive at the homogenious equation for a Linear Time-Invariant (LTI) dynamical system describing the pertubation, reading 

$$
\frac{\mathrm{d}\mathbf{q}'}{\mathrm{d}t}=\mathcal{L}(\mathbf{q}')
$$

with  

$$
\mathcal{L}\equiv\nabla_{\mathbf{q}}\mathcal{N}\big|_{\mathbf{q}_b}\in \mathbb{R}^{N\times N}
$$
representing the Jakobian evaluate at the base state  $\mathbf{q}_b$. 

### Spectral analysis of the  linear operator

We assume the pertubation to have the form of normal modes, reading
$$
\mathbf{q}'=\hat{\mathrm{q}}\mathrm{e^{-j\omega}t}+c.c.
$$
with the complex frequency $\omega = \omega_r+j\omega_i$.
Insertin thia in the linearized pertubation eqaution leads to ein eigenvalue problem, reading in matrix form 

$$
\mathbf{A} \hat{\mathbf{q}} = j \omega \hat{\mathbf{q}}
$$

where $A$ is the descrete version of the Jacobian $\mathcal{L}$,  $\omega$ is the eigenvalue and $\hat{q}$ is the eigenvector.

By solving the eigenvalue problem the imaginary part of the eigenvalues detmerines  

- $\omega_i>0$ := exponential temporal growth

- $\omega_i<0$ := exponential temporal decay

- $\omega_r$ := mode frequency
- $\hat{\mathbf{q}}$ :=  mode shape  


## FELiCS implementation (to be done)

The discretized eigenvalue equation is solved with a SLEPc-based solver in FELiCS. A set of eigenvalues $\omega_i$ ($i=1,2,...,n$) are calculated near a defined guess value $\sigma$.
For each $\omega_i$, the corresponding eigenvector $\hat{q_i}$ is also calculated.
The maximum number of $\omega_i$ is defined as $n$.
The real and imagine parts of the solved eigenvalues are gains and frequencies of specific mode while the real parts of the eigenvectors reveal the oscillation mode shape (or speak modal pattern).

### Derivation: For incompressible Navier-Stokes equations
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
$$
