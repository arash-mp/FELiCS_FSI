# Theory, governing equations and models

We introduce the math behind FELiCS and motivate the different types of linear analyses and the associated equations and models. The idea is that the reader understands why and for what cases the respective analyses, equations and models are used. Furthermore, the reader is supposed to get a brief overview of both the fundamental equations on a theory level and the actually implemented equations on a code level.

##### 1. Linear analysis types

We explain the specific linear analyses that the code is capable of solving, detailing its application scope. The theoretical background and motivation behind these analyses, which underpin the code’s implementation and use, are outlined.

```{toctree}
:maxdepth: 1
analysisMode_modal.md
analysisMode_resolvent.md
analysisMode_input_output.md
```

##### 2. Governing equations

We present the governing equations that describe the physical phenomena addressed by FELiCS. The nonlinear and linearized forms of these equations are detailed. The main assumptions and justifications used in deriving and formulating the equations are stated. Additionally, the motivation behind each governing equations is explained, and relevant use cases are illustrated with examples and references.

```{toctree}
:maxdepth: 1
equation_navier_stokes.md
equation_energy.md
equation_species.md
```

##### 3. Physical models and state equations

We describe the physical models and state equations required to close the governing equations in certain cases. The nonlinear and linearized forms of these models and equations are presented. The main assumptions and justifications underlying their derivation and formulation are outlined. Furthermore, the motivation for each model or equation is briefly sketched, accompanied by relevant use cases and examples.

```{toctree}
:maxdepth: 1
equation_state.md
equation_viscosity.md
equation_turbulence_models.md
equation_reaction.md
```

##### 4. Miscellaneous

We elaborate on additional components currently implemented in FELiCS, such as the sponge layer, providing an overview of their functionality and integration.

```{toctree}
:maxdepth: 1
equation_sponge.md
```