# Modal analysis
### This is the documentation for modal analysis operators. Xiuyang is modifying based on Thomas' and Sophie's comments. Welcome for any new reviewing and comments.

## Modal analysis for a Linear System
Consider the eigenproblem for a linear system:
$$
A \hat{q} = j \omega \hat{q}
$$
where $A$ is the linear operator, $j$ is the imaginary unit, $\omega$ is the eigenvalue and $\hat{q}$ is the eigenvector.

Applying the discretization scheme in a finite element method (FEM), a weighting metrix is added on the right hand side:
$$
A \hat{q} = \omega B \hat{q}
$$
where $B$ is the weighting matrix (with $j$ multiplied in).

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

## Residual calculation
After calculation, the relative residuals for each pair of eigenvalues and eigenvectors are calculated as
$$
\mathcal{R}_i = \frac{||A\hat{q_i} - \omega_i B \hat{q_i}||}{||\omega_i B \hat{q_i}||} .
$$
