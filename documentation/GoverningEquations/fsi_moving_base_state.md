# Rigid-body FSI about a mean-flow base state (`FSI.MovingBaseState`)

The rigid-body FSI coupling (`FELiCS/Equation/FSI/RigidBodyFSI.py`, after [Negi et al. 2020](#negi)) linearises the fluid and structure equations about a base state that is assumed to solve the steady equations FELiCS linearises. A time-averaged turbulent flow (DNS/LES mean) with a frozen eddy viscosity does not. With the option `FSI.MovingBaseState` the part of the mean momentum balance that the linear operator does not represent, and the eddy viscosity, move rigidly with the body. Then a body that is displaced together with its mean flow is again a steady solution, which the original formulation does not guarantee.

Assumptions:
- incompressible, $\rho = 1$, 2-D Cartesian geometry (same as the rest of the FSI coupling)
- the mean-flow residual and the eddy viscosity are attached to the body and follow its rigid motion
- the far field is not moved (inlet and lateral boundary conditions stay as they are)

Nomenclature:
- $\mathbf{U}_0, P_0$: base flow; $\mathbf{S}_0 = \tfrac12(\nabla\mathbf{U}_0 + \nabla\mathbf{U}_0^T)$
- $\mu = \mu_{lam} + \mu_{turb}$: dynamic viscosity used by FELiCS (`mu_tot`)
- $\eta_j, \phi_j = \dot\eta_j$: displacement and velocity of rigid-body DOF $j$
- $\mathbf{d}_j(\mathbf{x})$: displacement field of DOF $j$: $\mathbf{e}_j$ for a translation, $\hat{\mathbf{z}}\times(\mathbf{x}-\mathbf{c}_j)$ for a rotation about $\mathbf{c}_j$; $r_j = 1$ for a rotation and $0$ for a translation
- $\mathbf{n}$: outward normal of the fluid domain (it points into the body on the FSI interface $\Gamma$)
- $\mathbf{X}_u$: velocity test function

## Base-flow residual

FELiCS linearises the momentum equation with the frozen viscosity $\mu$,

$$
\mathcal{N}(\mathbf{U},P) = (\mathbf{U}\cdot\nabla)\mathbf{U} + \nabla P - \nabla\cdot(2\mu\mathbf{S}(\mathbf{U})) .
$$

A mean flow leaves the residual

$$
\mathbf{R}_0 = \mathcal{N}(\mathbf{U}_0,P_0) = \mathbf{C} - \mu_{lam}\mathbf{V} - \mathbf{T},\qquad
\mathbf{C} = (\mathbf{U}_0\cdot\nabla)\mathbf{U}_0 + \nabla P_0,\quad
\mathbf{V} = \nabla\cdot(2\mathbf{S}_0),\quad
\mathbf{T} = \nabla\cdot(2\mu_{turb}\mathbf{S}_0),
$$

i.e. the part of the Reynolds-stress force $-\nabla\cdot\overline{\mathbf{u}'\mathbf{u}'}$ that the eddy viscosity does not model. The mean flow is a steady solution of $\mathcal{N}(\mathbf{U},P) = \mathbf{R}_0(\mathbf{x})$. In the original formulation $\mathbf{R}_0$ is a force fixed in the laboratory frame. A body displaced through it feels a force proportional to the displacement, i.e. a spurious aerodynamic stiffness: a static heave of a body in a uniform stream produces lift.

## Moving base state

Let $\mathbf{R}_0$ and $\mu$ move rigidly with the body. For a rotation by $\eta_j$ a vector field transforms as $\mathbf{R}\mapsto\mathbf{Q}\mathbf{R}(\mathbf{Q}^T\mathbf{x})$, so to first order

$$
\delta\mathbf{R} = \eta_j\left[-(\mathbf{d}_j\cdot\nabla)\mathbf{R}_0 + r_j\,\hat{\mathbf{z}}\times\mathbf{R}_0\right],\qquad
\delta\mu = -\eta_j\,(\mathbf{d}_j\cdot\nabla)\mu .
$$

Linearising $\mathcal{N}_{\mu+\delta\mu}(\mathbf{U}_0+\mathbf{u}', P_0+p') = \mathbf{R}_0 + \delta\mathbf{R}$ gives the linearised momentum equation of FELiCS, $\mathcal{L}\mathbf{u}'$, with an additional forcing on the right-hand side:

$$
\partial_t\mathbf{u}' + \mathcal{L}\mathbf{u}' = \sum_j \eta_j\,\mathbf{f}_j,\qquad
\mathbf{f}_j = -(\mathbf{d}_j\cdot\nabla)\mathbf{R}_0 + r_j\,\hat{\mathbf{z}}\times\mathbf{R}_0 - \nabla\cdot\big(2\,m_j\,\mathbf{S}_0\big),\qquad
m_j = (\mathbf{d}_j\cdot\nabla)\mu .
$$

The force on DOF $k$ changes in two ways. The moving eddy viscosity changes the viscous traction by $2\,\delta\mu\,\mathbf{S}_0\mathbf{n}$. The base stress $\boldsymbol{\Sigma}_0 = P_0\mathbf{I} - 2\mu\mathbf{S}_0$ is evaluated on the displaced wall, which by the transport theorem for a surface moving with normal velocity $\mathbf{d}_j\cdot\mathbf{n}$ adds $(\mathbf{d}_j\cdot\mathbf{n})\,\mathbf{d}_k\cdot\nabla\cdot\boldsymbol{\Sigma}_0$, and $\nabla\cdot\boldsymbol{\Sigma}_0 = \mathbf{R}_0$ on the no-slip wall. Together

$$
F'_k = \int_\Gamma \mathbf{d}_k\cdot\big(p'\mathbf{n} - \boldsymbol{\tau}'\mathbf{n}\big)\,dS + \sum_j C_{kj}\,\eta_j,\qquad
C_{kj} = \underbrace{\int_\Gamma 2\,m_j\,\mathbf{d}_k\cdot(\mathbf{S}_0\mathbf{n})\,dS}_{G_{kj}}
+ \underbrace{\int_\Gamma (\mathbf{d}_j\cdot\mathbf{n})(\mathbf{d}_k\cdot\mathbf{R}_0)\,dS}_{Geo_{kj}} .
$$

$Geo_{kj}$ is the base-stress ("added stiffness") contribution discussed by [Negi et al. 2020](#negi), following Fanion et al. (2000). It vanishes at first order for a base flow that satisfies the steady equations at the wall ($\mathbf{R}_0 = 0$), which is why the original formulation omits it. For $\mathbf{R}_0 = 0$ and uniform $\mu$, all new terms vanish and the original formulation is recovered.

Check: for a static displacement, $\mathbf{u}' = -\eta_j[(\mathbf{d}_j\cdot\nabla)\mathbf{U}_0 - r_j\hat{\mathbf{z}}\times\mathbf{U}_0]$, $p' = -\eta_j(\mathbf{d}_j\cdot\nabla)P_0$ satisfies the forced equation exactly for any $\mathbf{U}_0$, $P_0$, $\mu$, and the transpiration condition on $\Gamma$. With it, $F'_k$ equals the force change of the rigidly displaced configuration: zero for a translation, and the rotation of the force vector for a rotation. Without the new terms this holds only if $\mathbf{R}_0 = 0$ and $\nabla\mu = 0$.

## Rotations and the wake (`MovingBaseStateWake`)

Rigid rotation is the right transport for the body and its boundary layer, but not for the wake. A pitched airfoil sheds its wake from the displaced trailing edge, and the wake then follows the free stream; it does not swing about the pivot. Rotating the residual rigidly by $\eta_j$ about $\mathbf{c}_j$ moves the wake at a distance $\ell$ behind the pivot by $\ell\,\eta_j$, which is several times the trailing-edge displacement a few chords downstream.

With `MovingBaseStateWake`, rotation DOFs carry $\mathbf{R}_0$ and $\mu$ with the field

$$
\boldsymbol{\xi}_j = w\,\mathbf{d}_j + (1-w)\,\mathbf{d}_j(\mathbf{x}_{TE}),\qquad
w = 1 - \mathrm{smoothstep}\!\left(\frac{x - x_{TE}}{L}\right),
$$

so the fields rotate rigidly for $x \le x_{TE}$ and are translated by the trailing-edge displacement for $x \ge x_{TE} + L$, with $L$ = `MovingBaseStateWakeLength`. The trailing edge $\mathbf{x}_{TE}$ is taken as the interface facet with the largest $x$. A vector field is carried with the local rotation of $\boldsymbol{\xi}_j$, $a_j = \tfrac12(\nabla\times\boldsymbol{\xi}_j)_z$, and without stretching:

$$
\mathbf{f}_j = -(\boldsymbol{\xi}_j\cdot\nabla)\mathbf{R}_0 + a_j\,\hat{\mathbf{z}}\times\mathbf{R}_0 - \nabla\cdot\big(2\,(\boldsymbol{\xi}_j\cdot\nabla\mu)\,\mathbf{S}_0\big).
$$

For a rigid rotation, $\boldsymbol{\xi}_j = \mathbf{d}_j$ and $a_j = 1$, which recovers the expression above. On the body $\boldsymbol{\xi}_j = \mathbf{d}_j$, so the force rows $G$ and $Geo$ are unchanged. Translations are not affected by the option.

```{note}
This is a modelling choice, not a derivation: the transported base state is no longer an exact steady solution in the blend region. The blend length $L$ is a parameter; check the sensitivity of the results to it.
```

## Weak form and discrete coupling

The forcing is integrated by parts, with $\boldsymbol{\xi}_j = \mathbf{d}_j$, $a_j = r_j$ and $\nabla\cdot\boldsymbol{\xi}_j = 0$ unless `MovingBaseStateWake` is on. As in the momentum equation, the viscous boundary term is dropped:

$$
\begin{aligned}
g_j(\mathbf{X}_u) = {} & \int_\Omega \mathbf{R}_0\cdot\big((\boldsymbol{\xi}_j\cdot\nabla)\overline{\mathbf{X}_u}\big)\,dx
+ \int_\Omega (\nabla\cdot\boldsymbol{\xi}_j)\,\mathbf{R}_0\cdot\overline{\mathbf{X}_u}\,dx
- \oint_{\partial\Omega}(\boldsymbol{\xi}_j\cdot\mathbf{n})\,\mathbf{R}_0\cdot\overline{\mathbf{X}_u}\,dS \\
& + \int_\Omega a_j\,(\hat{\mathbf{z}}\times\mathbf{R}_0)\cdot\overline{\mathbf{X}_u}\,dx
+ \int_\Omega 2\,(\boldsymbol{\xi}_j\cdot\nabla\mu)\,\mathbf{S}_0 : \nabla\overline{\mathbf{X}_u}\,dx .
\end{aligned}
$$

With FELiCS' convention $A = -\mathrm{j}\times$(algebraic terms), the augmented operator gets

- $A[\text{momentum rows}, \eta_j] = +\mathrm{j}\,g_j$, except in Dirichlet rows and in the interface rows that hold the transpiration constraint;
- $A[\phi_k, \eta_j] = -\mathrm{j}K_{kj} + \mathrm{j}C_{kj}$.

$B$ is unchanged.

## Input

`MovingBaseStateResidual` selects how $\mathbf{R}_0$ is obtained.

- `"Discrete"` (default): FELiCS computes the residual of its own P2 base flow under its own momentum operator,

  $$
  (\mathbf{R}_0, \mathbf{v}) = \int_\Omega \big[(\nabla\mathbf{U}_0)\,\rho\mathbf{U}_0\big]\cdot\mathbf{v} + \boldsymbol{\tau}_0 : \nabla\mathbf{v} - P_0\,\nabla\cdot\mathbf{v}\,dx
  + \oint_{\partial\Omega} (P_0\mathbf{n} - \boldsymbol{\tau}_0\mathbf{n})\cdot\mathbf{v}\,dS
  \quad \forall\,\mathbf{v}\in P_2,
  $$

  with $\boldsymbol{\tau}_0 = \mu(\nabla\mathbf{U}_0 + \nabla\mathbf{U}_0^T)$. This is the L2 projection of the strong residual including the inter-element flux jumps, i.e. exactly the force the discrete linear operator treats as balancing the base flow. No extra input is needed.
- `"File"`: the pieces $\mathbf{C}$, $\mathbf{V}$, $\mathbf{T}$ are read from the mean-flow file as the fields `fsiCx, fsiCy, fsiVx, fsiVy, fsiTx, fsiTy`, computed from the source data (e.g. spectral derivatives of DNS statistics) with the same $\mu_{turb}$ as `muturb`/`nuturb`. FELiCS assembles $\mathbf{R}_0 = \mathbf{C} - \mu_{lam}\mathbf{V} - \mathbf{T}$, where $\mathbf{T}$ is included only for `TurbulenceModel = "File"`. A missing field stops the run with an error.

```{warning}
In a boundary layer that the P2 mesh does not fully resolve, the residual of the continuous mean flow is the small difference of large convective and viscous terms, and it can differ by O(1) from what the discrete operator sees. Moving the "File" residual then leaves part of the spurious stiffness (for a NACA0017 mean flow at $Re = 1.6\times10^5$, a static heave force of +2.3 per unit displacement instead of about 0, versus −5.7 without the option). The "Discrete" residual is consistent with the operator by construction and gave −0.2.
```

```json
"FSI": {
    "Enabled": true,
    "MovingBaseState": true,
    "MovingBaseStateResidual": "Discrete",
    "StaticCheck": true,
    ...
}
```

- `MovingBaseState` (default `false`): switches the terms above on.
- `MovingBaseStateResidual` (default `"Discrete"`): source of $\mathbf{R}_0$, `"Discrete"` or `"File"`.
- `MovingBaseStateWake` (default `false`, needs `MovingBaseState`): rotation DOFs use the wake-following transport instead of a rigid rotation everywhere.
- `MovingBaseStateWakeLength` (default `0.5`, chord units): blend length $L$ downstream of the trailing edge.
- `MovingBaseStateForceTerms` (default `"All"`, needs `MovingBaseState`): which force-row terms are added, for sensitivity tests. `"All"` = $G + G^{eo}$; `"NoGeo"` drops the base stress at the displaced wall $G^{eo}$ (the "added stiffness" that Negi et al. 2020 show to vanish for an exact base flow); `"NoG"` drops the moving eddy-viscosity traction $G$; `"None"` drops both. The momentum forcing is always kept. Both matrices are written to the log in every case. Note that only `"All"` keeps the exact result that a static translation gives zero force.
- `StaticCheck` (default `false`): solves the steady fluid problem for a unit static displacement of every DOF and logs the aerodynamic force matrix, for the frozen base state and (if enabled) the moving one. A pure translation in a uniform stream should give nearly zero. This costs one extra LU factorisation.

## Per-mode output

For every FSI mode the log reports the force on each DOF and its split into an apparent aerodynamic stiffness and damping. With $s = -\mathrm{j}\lambda$ (perturbations $\propto e^{-\mathrm{j}\lambda t}$):

$$
\frac{F'_k}{\eta_k} = -k_{a,k} - c_{a,k}\,s,\qquad k_{a,k},\,c_{a,k}\in\mathbb{R}.
$$

For a single DOF this gives $k_a = m|\lambda|^2 - K$ and $c_a = -2m\,\mathrm{Im}(\lambda) - D$; a mode grows when $c_a < -D$. For several DOFs, $F'_k/\eta_k$ includes the cross-coupling, so $k_a$ and $c_a$ are apparent per-DOF values. The log also gives the relative residual of the structural equation, which should be at solver accuracy. Modes with a structural energy ratio $\geq 0.1$ are logged at info level, the rest at debug level.

### References

<a id="negi">[1]</a> P. S. Negi, A. Hanifi, D. S. Henningson, "On the linear global stability analysis of rigid-body motion fluid–structure-interaction problems", J. Fluid Mech. 903, A35 (2020). Preprint: [arXiv:1910.09605](https://arxiv.org/abs/1910.09605)
