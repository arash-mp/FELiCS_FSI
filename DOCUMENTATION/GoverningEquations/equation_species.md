# Species transport equations
The mass/ species transport equations describe the transport behaviours of specific component of the mixtures, or more general, a passive scalar field.

Motivation:
- Component concentration can affect phisical properties, so it needs to be track when a mixture is applied;
- The transport of passive scalar fields are very important in some cases. For example progress variables in combustion models.

Use case:
- in flows with chemical reactions;
- in flows with mixture which has non-uniform concentration distribution;
- tutorial example: passive scalar transport, turbulent flame?

References:
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)
- [Kaiser et al. 2020](https://doi.org/10.1017/jfm.2021.151)


Assumptions:
- The density $\rho$ usually changes due to chemical reactions, changes of component concentrations, high-Mach numbers, etc...
- Source terms can be defined based on reaction / combustion models
- no gravitational forces

### Nonlinear equations
$$
\rho \frac{\partial Y}{\partial t} + \rho \mathbf{u} \cdot \nabla Y = \nabla \cdot (\rho D\nabla Y) + f_{Y}
$$
where $D$ is mass diffusivity, and $f_Y$ is the source term.

__The nonlinear equations are currently not implemented in FELiCS.__



### Linear equations
$$
\overline{\rho} \frac{\partial Y'}{\partial t} + \rho ' \overline{\mathbf{u}} \cdot \nabla \overline{Y} + \overline{\rho} \mathbf{u}' \cdot \nabla \overline{Y} + \overline{\rho} \overline{\mathbf{u}} \cdot \nabla Y' + \rho'\mathbf{u}' \cdot \nabla Y' = \nabla \cdot \rho' \overline{D} \nabla \overline{Y} + \nabla \cdot \overline{\rho} D' \nabla \overline{Y} + \nabla \cdot \overline{\rho} \overline{D} \nabla Y' + f_{Y}'
$$
where the last term on the left hand side $\rho' \mathbf{u} \cdot \nabla Y'$ is the nonlinear term and need to be modeled with additional equations within turbulent flows.

#### Weak form
The weak form of the linearized species transport equation with normal mode ansatz, as implemented in FELiCS, is
$$
\int \omega \overline{\rho} \widehat{Y}_s \bf{X}^* dx =\int j \, \widehat{Y}_s \, \nabla \cdot \left( \overline{\rho} \, \overline{\bf{u}} \, \bf{X}^* \right) \, dx+ \int j \, \overline{Y}_s \, \nabla \cdot \left( \overline{\rho} \, \widehat{\bf{u}} \, \bf{X}^* \right) \, dx+ \int j \, \overline{Y}_s \, \nabla \cdot \left( \widehat{\rho} \, \overline{\bf{u}} \, \bf{X}^* \right) \, dx - \int j \, \widehat{Y}_s \, \overline{\rho} \, \overline{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds- \int j \, \overline{Y}_s \, \overline{\rho} \, \widehat{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds- \int j \, \overline{Y}_s \, \widehat{\rho} \, \overline{\bf{u}} \cdot \bf{X}^* \cdot \bf{n} \, ds - \int j \, \overline{D}_s \, \nabla \widehat{Y}_s \cdot \nabla \bf{X}^* \, dx- \int j \, \widehat{D}_s \, \nabla \overline{Y}_s \cdot \nabla \bf{X}^* \, dx + \int \overline{f}_s \, \bf{X}^* \, dx
$$

