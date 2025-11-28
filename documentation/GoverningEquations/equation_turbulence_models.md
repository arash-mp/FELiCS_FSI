# Turbulence models

Turbulence models in FELiCS are used to represent the momentum transfer caused by unresolved turbulent motions through an eddy viscosity added to the molecular viscosity. In the mean (Reynolds- or Favre-averaged) equations this models the deviatoric part of the Reynolds stresses. In the linearized fluctuation equations it augments the diffusive action on perturbations.

```{note}
In the current FELiCS implementation, the eddy viscosity is frozen: only a prescribed mean field $\overline{\mu_t(\mathbf{x})}$ (or a spatially constant) is used, and no fluctuation of eddy viscosity is modeled.
Linearized turbulence models are scheduled for a later release.
```

Assumptions:
- turbulent Reynolds stresses are proportional to the mean strain rate
- turbulent mixing is analogous to molecular viscosity, using an effective eddy viscosity
- turbulent viscosity is isotropic

Example use case:
- resolvent analysis of turbulent jets: [von Saldern et al. 2024](https://doi.org/10.1017/jfm.2024.922)
- turbulent swirling jet: [Müller et al. 2020](https://doi.org/10.1017/jfm.2019.1063)
- turbulent flow in a stenosis: [Resolvent analysis](./../Tutorials/Resolvent_Analysis.md)

Nomenclature:
- $\mathbf{u}$: velocity  
- $p$: pressure  
- $\rho$: density  
- $\overline{(\cdot)}$, $(\cdot)'$: mean and fluctuation (Reynolds or Favre average, as appropriate)  
- $\boldsymbol{\tau}$: viscous stress tensor  
- $\mu$: molecular (dynamic) viscosity
- $\mu_t$: eddy (turbulent) viscosity  
- $\mu_\text{eff} = \mu + \mu_t$: effective viscosity  
- $k$: turbulent kinetic energy
- $\mathbf{I}$: identity tensor  
- $\overline{\mathbf{S}}=\tfrac12(\nabla\overline{\mathbf{u}}+\nabla\overline{\mathbf{u}}^{T})$, $\mathbf{S}'=\tfrac12(\nabla\mathbf{u}'+\nabla\mathbf{u}'^{T})$: rate of strain tensor
- $\rho \overline{\mathbf{u}'\mathbf{u}'}$: Reynolds stress tensor

References:
- [Boussinesq 1872](https://gallica.bnf.fr/ark:/12148/bpt6k56673076/f2.item.texteImage)
- [Pope 2000](https://doi.org/10.1017/CBO9781316179475)
- [Reynolds & Hussain 1972](https://doi.org/10.1017/S0022112072000679)

## Theoretical background

The following illustrates the use of an eddy viscosity in ther Navier--Stokes equations for an incompressible flow.

### Mean equations and mean eddy viscosity
After applying the Reynolds decomposition to the velocity and pressure, the mean momentum equation reads

$$
\rho(\overline{\mathbf{u}}\!\cdot\!\nabla)\overline{\mathbf{u}}
= -\nabla\overline{p} + \mu\nabla^2\overline{\mathbf{u}}
- \rho\,\nabla\!\cdot\overline{\mathbf{u}'\mathbf{u}'}.
$$

Under the Boussinesq hypothesis, the mean Reynolds stress is modeled as

$$
-\rho\,\overline{\mathbf{u}'\mathbf{u}'} + \tfrac{2}{3}\rho k\,\mathbf{I}
\;\approx\; 2\,\mu_t\,\overline{\mathbf{S}}.
$$

so that the deviatoric part of $-\rho\,\overline{\mathbf{u}'\mathbf{u}'}$ is represented by an eddy viscosity $\mu_t=\rho\nu_t$. The isotropic part $-\tfrac{2}{3}\rho k\,\mathbf{I}$ can be absorbed into a modified mean pressure. Hence, the mean viscous stress becomes

$$
\overline{\boldsymbol{\tau}} \;=\; \overline{\mu}_\text{eff}\,(\nabla+\nabla^T)\overline{\mathbf{u}},
\qquad
\overline{\mu}_\text{eff} \equiv \mu + \mu_t(\mathbf{x}).
$$

### Fluctuation equations and fluctuating Reynolds stress

Subtracting the mean from the instantaneous equations gives, for incompressible flows,

$$
\rho\Big(\partial_t \mathbf{u}'
+ \overline{\mathbf{u}}\!\cdot\!\nabla \mathbf{u}'
+ \mathbf{u}'\!\cdot\!\nabla \overline{\mathbf{u}}\Big)
= -\nabla p' + \mu\nabla^2\mathbf{u}' \;-\; \rho\,\nabla\!\cdot\Big[\mathbf{u}'\mathbf{u}'-\overline{\mathbf{u}'\mathbf{u}'}\Big].
$$

We define the fluctuating Reynolds stress,

$$
\mathbf{r}' \equiv \mathbf{u}'\mathbf{u}' - \overline{\mathbf{u}'\mathbf{u}'}.
$$

A linearized Boussinesq closure reads schematically

$$
\boldsymbol{\tau}'_{\text{turb}} \;\approx\; 2\,\overline{\mu}_t\,\mathbf{S}' \;+\; 2\,\mu_t'\,\overline{\mathbf{S}} \;-\; \tfrac{2}{3}\rho k'\,\mathbf{I}.
$$

where
- the $2\,\overline{\mu}_t\,\mathbf{S}'$ term models turbulent diffusion due to velocity fluctuations;
- the $2\,\mu_t'\,\overline{\mathbf{S}}$ term models turbulent diffusion due to eddy viscosity fluctuations;
- the isotropic parts can again be absorbed into $p'$.

```{note}
If compressibility is considered, the stress form used in FELiCS mirrors the molecular one$\boldsymbol{\tau}_\text{eff} = \mu_\text{eff}\big[(\nabla+\nabla^T)\mathbf{u} - \tfrac{2}{3}(\nabla\!\cdot\!\mathbf{u})\mathbf{I}\big]$, with $\mu_\text{eff}=\mu+\mu_t$.
```

```{note}
Currently, the fluctuations of eddy viscosity are neglected $(\mu_t'=0)$ and not implemented in FELiCS.
```

```{note}
There is currently no transport/closure equation for $\overline{\mu}_t$ being solved in FELiCS.
```


## Turbulence models currently implemented in FELiCS

### Constant eddy viscosity $(\overline{\mu}_t= \textrm{const})$

A single scalar applied uniformly in space. To use this option, directly set an effective viscosity value $\overline{\mu}_\textrm{eff}=\overline{\mu}_t + \overline{\mu}$ into the molecular viscosity field of the FELiCS [configuration file](./../Running_FELiCS/FELiCS_settings.md):

```yaml
  "MolViscModel": "Constant"
  "MolVisc": (float)
```

For a laminar flow with a constant viscosity, simply set the molecular value to the `"MolVisc"` entry. 

### Arbitrary eddy viscosity scalar field from file $(\overline{\mu}_t(\mathbf{x}))$

Alternatively, an arbitrary field defined on the import mesh can be considered. This can be useful, for example, when importing the mean flow eddy viscosity from an outside RANS solver. 

In such case, the FELiCS [configuration file](./../Running_FELiCS/FELiCS_settings.md) should contain:

```yaml
  "TurbulenceModel": "File"
```

and the [mean flow file](./../Running_FELiCS/fel_file.md) should contain an array called `nuturb`, defining $\overline{\mu}_t(\mathbf{x})$ over the mean flow mesh.

```{note}
Small FELiCS quirk: Even though the eddy and molecular viscosity are called `nuturb` and `nulam` in the code, they actually refer to the dynamic viscosity $\mu$.
```

```{note}
Currently, the $k$-$\varepsilon$ equations are being implemented and validated in FELiCS and will soon be publicly released.
```