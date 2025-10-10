# Turbulence models

Turbulence models in FELiCS are used to represent the momentum transfer caused by unresolved turbulent motions through an eddy viscosity added to the molecular viscosity. In the mean (RANS/Favre-averaged) equations this models the deviatoric part of the Reynolds stress; in the linearized fluctuation equations it augments the diffusive action on perturbations.

```{note}
In the current FELiCS implementation, the eddy viscosity is **frozen**: only a prescribed mean field $\overline{\mu_t(x)}$ (or a constant) is used, and no fluctuation of eddy viscosity is modeled.
Linearized turbulence models are scheduled for a later release.
```

Assumptions:
- turbulent Reynolds stresses are proportional to the mean strain rate
- turbulent mixing is analogous to molecular viscosity, using an effective eddy viscosity
- turbulent viscosity is isotropic

Example use case:
- Resolvent analysis of turbulent jets: [von Saldern et al. 2024](https://doi.org/10.1017/jfm.2024.922)
- Turbulent swirling jet: [Müller et al. 2020](https://doi.org/10.1017/jfm.2019.1063)
- Turbulent flow in a stenosis: [Resolvent analysis](./../Tutorials/Resolvent_Analysis.md)

Nomenclature:
- $\mathbf{u}, p, \rho$: velocity, pressure, density  
- $\overline{(\cdot)}$, $(\cdot)'$: mean and fluctuation (Reynolds or Favre as appropriate)  
- $\boldsymbol{\tau}$: viscous stress tensor  
- $\mu$: molecular (dynamic) viscosity; 
- $\mu_t$: eddy (turbulent) viscosity  
- $\mu_\text{eff} = \mu + \mu_t$: effective viscosity  
- $k$: turbulent kinetic energy; 
- $\mathbf{I}$: identity tensor  
- $\overline{\mathbf{S}}=\tfrac12(\nabla\overline{\mathbf{u}}+\nabla\overline{\mathbf{u}}^{T})$, $\mathbf{S}'=\tfrac12(\nabla\mathbf{u}'+\nabla\mathbf{u}'^{T})$: Rate of strain tensor
- $\rho \overline{\mathbf{u}'\mathbf{u}'}$: Reynolds stress tensor

References:
- [Boussinesq 1872](https://gallica.bnf.fr/ark:/12148/bpt6k56673076/f2.item.texteImage)
- Pope (2000), Turbulent Flows
- Wilcox (2006), Turbulence Modeling for CFD
- Reynolds & Hussain (1972), J. Fluid Mech.

## Theoretical background

The following illustrates the use of an eddy viscosity in ther Navier-Stokes equations for an incompressible flow.

### Mean equations and mean eddy viscosity
After applying the Reynold's decomposition to the velocity and pressure, the mean momentum equation reads

$$
\rho(\overline{\mathbf{u}}\!\cdot\!\nabla)\overline{\mathbf{u}}
= -\nabla\overline{p} + \mu\nabla^2\overline{\mathbf{u}}
- \rho\,\nabla\!\cdot\overline{\mathbf{u}'\mathbf{u}'}.
$$

Under the **Boussinesq hypothesis**, the mean Reynolds stress is modeled as

$$
-\rho\,\overline{u_i' u_j'} + \tfrac{2}{3}\rho k\,\delta_{ij}
\;\approx\; 2\,\mu_t\,\overline{S}_{ij}.
$$

so that the **deviatoric** part of $-\rho\,\overline{\mathbf{u}'\mathbf{u}'}$ is represented by an **eddy viscosity** $\mu_t=\rho\nu_t$. The isotropic part $-\tfrac{2}{3}\rho k\,\mathbf{I}$ is absorbed into a modified mean pressure.   Hence, the mean viscous stress becomes

$$
\overline{\boldsymbol{\tau}} \;=\; \overline{\mu}_\text{eff}\,(\nabla+\nabla^T)\overline{\mathbf{u}},
\qquad
\overline{\mu}_\text{eff} \equiv \mu + \mu_t(\mathbf{x}).
$$

### Fluctuation equations and fluctuating Reynolds stress

Subtracting the mean from the instantaneous equations gives (incompressible)

$$
\rho\Big(\partial_t \mathbf{u}'
+ \overline{\mathbf{u}}\!\cdot\!\nabla \mathbf{u}'
+ \mathbf{u}'\!\cdot\!\nabla \overline{\mathbf{u}}\Big)
= -\nabla p' + \mu\nabla^2\mathbf{u}' \;-\; \rho\,\nabla\!\cdot\Big[\mathbf{u}'\mathbf{u}'-\overline{\mathbf{u}'\mathbf{u}'}\Big].
$$

Define the **fluctuating Reynolds stress**,

$$
\mathbf{r}' \equiv \mathbf{u}'\mathbf{u}' - \overline{\mathbf{u}'\mathbf{u}'}.
$$

A linearized, viscosity-form Boussinesq closure (written directly for **stress**) reads schematically

$$
\boldsymbol{\tau}'_{\text{turb}} \;\approx\; 2\,\overline{\mu}_t\,\mathbf{S}' \;+\; 2\,\mu_t'\,\overline{\mathbf{S}} \;-\; \tfrac{2}{3}\rho k'\,\mathbf{I}.
$$

where:
- The $2\,\overline{\mu}_t\,\mathbf{S}'$ term **augments diffusion of the fluctuations** (a viscous-like action on $\mathbf{u}'$ weighted by $\overline{\mu}_t$).  
- The $2\,\mu_t'\,\overline{\mathbf{S}}$ term would couple **fluctuations of eddy viscosity** to the **mean strain** (currently not implemented in FELiCS).  
- Isotropic parts can again be absorbed into $p'$.

When $\overline{\mu}_t$ is **uniform** (and also $\overline{\mu}$), the linear viscous contribution reduces to $(\overline{\mu}+\overline{\mu}_t)\nabla^2\mathbf{u}'$. For spatially varying $\overline{\mu}_t(\mathbf{x})$, the operator includes $\nabla\!\cdot\!\big(2\,\overline{\mu}_t\,\mathbf{S}'\big)$ with the appropriate product-rule terms.

```{note}
**Compressible note:** if compressibility is considered, the stress form used in FELiCS mirrors the molecular one:  
$$
\boldsymbol{\tau}_\text{eff} = \mu_\text{eff}\big[(\nabla+\nabla^T)\mathbf{u} - \tfrac{2}{3}(\nabla\!\cdot\!\mathbf{u})\mathbf{I}\big],
$$  
with $\mu_\text{eff}=\mu+\mu_t$.
```

## FELiCS-specific hypotheses
- **Frozen eddy viscosity:** Currently, the fluctuations of eddy viscosity are neglected $(\mu_t'=0)$.  
- **No mean eddy-viscosity model:** there is currently **no** transport/closure equation for $\overline{\mu}_t$ being solved in FELiCS.

## Turbulence models currently implemented in FELiCS

### Constant eddy viscosity $(\overline{\mu}_t=cst.)$

A single scalar applied uniformly. To use this option, directly set an effective viscosity value $\overline{\mu}_{eff}=\overline{\mu}_t + \overline{\mu}$ into the **molecular viscosity** field of the FELiCS [configuration file](./../Running_FELiCS/FELiCS_settings.md):

```yaml
  "MolViscModel": "Constant"
  "MolVisc": $mu_eff$,
```

For a laminar flow with a constant viscosity, simply set the molecular value to the `"MolVisc"` entry. 

### Arbitrary eddy viscosity scalar field from file $(\overline{\mu}_t(x))$

Alternatively, an arbitrary field defined on the import mesh can be considered. This can be useful, for example, when importing the mean flow eddy viscosity from an outside RANS solver. 

In such case, the FELiCS [configuration file](./../Running_FELiCS/FELiCS_settings.md) should contain:

```yaml
  "TurbulenceModel": "File"
```

and the [mean flow file](./../Running_FELiCS/fel_file.md) should contain an array called `nuturb`, defining $\overline{\mu}_t(x)$ over the mean flow mesh.

```{note}
:fearful: Small FELiCS quirck: eventhough the eddy and molecular viscosity are called `nuturb` and `nulam` in the code, they actually refer to the dynamic viscosity $\mu$. :fearful:
```

```{note}
Currently, the $k$-$\varepsilon$ equations are being implemented and validated in FELiCS and will soon be publicly released.
```