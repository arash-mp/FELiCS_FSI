# Theory and implemented equations
We introduce the math behind FELiCS and motivate the different types of linear analyses and the associated equations and models. The idea is that the reader understands why and for what cases the respective analyses, equations and models are used. Furthermore, the reader is supposed to get a brief overview of both the fundamental equations on a theory level and the actually implemented equations on a code level.

##### 1. Linear analysis types

- What can the code solve for?
- What is the theory and motivation behind these types of analyses?

```{toctree}
:maxdepth: 1
analysisMode_modal.md
analysisMode_resolvent.md
analysisMode_input_output.md
```

##### 2. Partial differential equations

- What are the governing equations for the physics?
- What are the nonlinear, linear (and bilinear, or second-order derivative) forms of the equations?
- What are the main assumptions and justifications used in deriving/formulating the equations?
- What is the motivation or what are the use cases for each equations (with some references/examples)?

```{toctree}
:maxdepth: 1
equation_navier_stokes.md
equation_energy.md
equation_species.md
```

##### 3. Algebraic equations

- These equations are required for closing the PDEs in some cases.
- Nonlinear and linear (and bilinear maybe) forms of the equations.
- What are the main assumptions and justifications used in deriving/formulating the equations?
- What is the motivation or what are the use cases for each equations (with some references/examples)?

```{toctree}
:maxdepth: 1
equation_state.md
equation_viscosity.md
equation_reaction.md
```

##### 4. Miscellaneous

- All other stuff (for now only sponge and tensor framework)

```{toctree}
:maxdepth: 1
equation_sponge.md
misc_tensor_formalism.md
```