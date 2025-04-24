# Implementation of modal analysis operators
### This is the documentation for modal analysis operators. Xiuyang has done the draft. Welcome for reviewing and comments.

## Modal analysis for a Linear System
Consider the eigenproblem for a linear system:
$$
A \hat{q} = j \omega \hat{q}
$$
where $A$ is the linear operator, $j$ is the imaginary unit, $\omega$ is the eigenvalue and $\hat{q}$ is the eigenvector.

Applying the discretization scheme in a finite element method (FEM), a weighting metrix is added on the right hand side.
$$
A \hat{q} = \omega B \hat{q}
$$
where $B$ is the weghting metrix (with $j$ multiplied in).

The discretized eigenvalue equation is solved with a SLEPc-based solver in FELiCS. A set of eigenvalues $\omega_i$ ($i=1,2,...,n$) can be get near a defined guess value $\sigma$.
For each $\omega_i$, the corrisponding eigenvector $\hat{q_i}$ is also calculated.
The maximum number of $\omega_i$ is defined as $n$.

## Residual calculation
After calculation, the relative residuals for each pair of eigenvalues and eigenvectors are calculated as
$$
\mathcal{R}_i = \frac{||A\hat{q_i} - \omega_i B \hat{q_i}||}{||\omega_i B \hat{q_i}||} 
$$
