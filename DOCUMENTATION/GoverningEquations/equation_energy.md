# Energy equations
The energy equations are usually required for compressible flows when thermodynamic changes (e.g. changes in enthalpy or temperature) in the flow field become relevant, such as flows with acoustics or chemical reactions. There are a lot of different forms of energy equations, depending on the given physical problem and the made associated assumptions.

Motivation:
- fundamental equations for flows with thermodynamic processes
- conservation of energy (or some form of it)

Use case:
- flows with acoustics, flows with combustion/chemical reactions
- tutorial example: turbulent flame? heat conduction?

References:
- Demange2024 AIAA
- Kaiser2023 CnF

Nomenclature:
- $\mathbf{u}$: velocity vector
- $p$: pressure
- $T$: temperature
- $T_b$: burnt, adiabatic flame temperature
- $T_u$: unburnt temperature
- $\rho$: density
- $\gamma$: heat capacity ratio
- $\kappa$: thermal conductivity
- $\tau$: viscous stress tensor
- $\mu$: dynamic viscosity
- $c$: Progress variable


## Energy pressure equation
Assumptions:
- primitive form
- constant Prandtl number
- constant heat capacity ratio
- molecular viscosity adheres to Sutherland's law
- How does thermal conductivity change?

### Nonlinear equations
$$
\frac{\partial p}{\partial t} + \mathbf{u} \cdot \nabla p + \gamma p(\nabla \cdot \mathbf{u}) - (\gamma - 1)\left[ \nabla (\kappa \nabla T) + \tau : \nabla \mathbf{u}\right] = 0
$$

with the mean viscous shear stress tensor $\overline{\tau}$ being defined as

$$
\overline{\tau} = \overline{\mu} \left[\left(\nabla + \nabla ^T\right)\overline{\mathbf{u}} - \frac{2}{3} (\nabla \cdot \overline{\mathbf{u}}) \mathbf{I}\right]
$$

### Linear equations
$$
\frac{\partial p'}{\partial t} + \overline{\mathbf{u}} \cdot \nabla p' + \mathbf{u}' \cdot \nabla \overline{p} + \gamma \left[ \overline{p}(\nabla \cdot \mathbf{u}') + p'(\nabla \cdot \overline{\mathbf{u}}) \right] - (\gamma - 1)\left[ \nabla (\overline{\kappa} \nabla T' + \kappa' \nabla \overline{T} + \overline{\tau} : \nabla \mathbf{u}' + \tau' : \nabla \overline{\mathbf{u}}) \right] = 0
$$

with the fluctuating viscous shear stress tensor $\tau'$ being defined as

$$
\tau' = \mu' \left[\left(\nabla + \nabla ^T\right)\overline{\mathbf{u}} - \frac{2}{3} (\nabla \cdot \overline{\mathbf{u}}) \mathbf{I}\right] + \overline{\mu}\left[(\nabla + \nabla ^T) \mathbf{u}' - \frac{2}{3} (\nabla \cdot \mathbf{u}') \mathbf{I}\right]
$$
with $\mathbf{I}$ being the identity tensor.





## Energy enthalpy equation
Assumptions:


### Nonlinear equations


### Linear equations




## Progress variable equation
Motivation and use case:
- active flame model (Kaiser CnF 2023) = fluctuations of reaction rate are included
- relate progress variable with reaction rate
- very simple but already effective active flame model under turbulent conditions

Assumptions/justifications:
- progress variable and density are directly related through algebraic equations, which obviates the use of an explicit energy equation $\rightarrow$ the species equation also acts as an energy conservation equation
- incompressible flow due to (1) low Mach number, (2) flame being acoustically compact
- adiabatic walls (no heat losses), adiabatic flame temperature is reached in the entire domain
- constant heat capacity and constant specific gas constants
- all assumptions result in temperature only being a function of the progress variable, $T = T_u + (T_b - T_u)c$
- low Mach equation of state, pressure is constant [link to equation of state page for low Mach]


### Nonlinear equations
[link to species qeuation page]


### Linear equations
[link to species qeuation page]