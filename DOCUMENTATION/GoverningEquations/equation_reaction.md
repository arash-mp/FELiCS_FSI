# Reaction models
The reaction models are used when chemical reactions due to a flame occur in the flow. The terms of the reaction models typically occur as source terms in the energy or species transport equation.


## Eddy break-up model
In the modeling of turbulent flames within a RANS framework, turbulence-flame interactions are linked through a characteristic turbulent time scale. The reaction rate of the flame is directly related to the state variables of the turbulent flow, allowing a coupled representation of combustion and turbulence dynamics.

Assumptions:
- turbulence controls combustion rate, reaction rate is proportional to the turbulence dissipation rate
- mixing-limited combustion, chemical reaction is assumed fast, mixing of reactants limits the combustion rate
- flame structure not resolved, no detailed flame structure is modeled, averaged quantities are used
- local equilibrium assumed, reactants are assumed to be locally in chemical equilibrium once mixed

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
where $A_\textrm{EBU}$ is a model constant lumped together with the turbulent time scale $\overline{\varepsilon}/\overline{k}$.


### Linearized eddy break-up model
The reaction rate of the linearized eddy break-up model is
$$
\dot{\Omega}' = A_{EBU} (\rho' (\overline{c} - \overline{c}^2) + \overline{\rho}(c'-2\overline{c}c'))
$$
with the same model constant $A_\textrm{EBU}$ as for the mean reaction rate.