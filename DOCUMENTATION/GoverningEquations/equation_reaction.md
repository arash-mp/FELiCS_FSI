# Reaction models


## Eddy break-up model

Assumptions:
- For the assumptions of the RANS EBU model the interested reader is kindly referred to [Poinsot & Veynante 2005](http://refhub.elsevier.com/S0010-2180(23)00162-1/sbref0058).

Example use case:
- Reacting turbulent jet flame: [Input-output analysis](./../Tutorials/input_output_analysis.md)

References:
- [Poinsot & Veynante 2005](http://refhub.elsevier.com/S0010-2180(23)00162-1/sbref0058)
- [Kaiser et al. 2023](https://doi.org/10.1016/j.combustflame.2023.112778)

Nomenclature:
- $\dot{\Omega}$: chemical reaction rate
- $C_\textrm{EBU}$: model constant
- $A_\textrm{EBU}$: model constant lumped together turbulent time scale $\overline{\varepsilon}/\overline{k}$
- $k$: turbulent kinetic energy
- $\varepsilon$: turbulent dissipation
- $\rho$: density
- $c$: progress variable

The mean reaction rate is
$$
\overline{\dot{\Omega}} = C_{EBU} \frac{\overline{\varepsilon}}{\overline{k}} \overline{\rho} \overline{c} (1-\overline{c}) =  A_{EBU} \overline{\rho} \overline{c} (1-\overline{c})
$$
where $A_\textrm{EBU}=C_{EBU} \frac{\overline{\varepsilon}}{\overline{k}}$ is constant in time, but may vary in space. It can be determined by inserting the temporal mean state, obtained by e.g. LES, into the equation above. 


### Linearized eddy break-up model
The reaction rate of the linearized eddy break-up model is
$$
\dot{\Omega}' = A_{EBU} (\rho' (\overline{c} - \overline{c}^2) + \overline{\rho}(c'-2\overline{c}c'))
$$
with the same model constant $A_\textrm{EBU}$ as for the mean reaction rate.
