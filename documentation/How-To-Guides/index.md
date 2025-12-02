# How-To-Guides

The FELiCS **How-To** guides are designed as an entry point for exploring **how to use FELiCS as a Python package**. Each guide is based on Jupyter Notebooks, making them interactive and easy to follow. These example also demonstrate how FELiCS can be integrated into Python projects.

All notebooks, together with the required files, are located in the [`how-to` folder]
(https://gitlab.com/felics-group/FELiCS/-/tree/main/how_tos?ref_type=heads).

Before working through the guides, **it is essential to install FELiCS as a package**. Detailed instructions for installation can be found in
the [installation guide](https://felics-d43476.gitlab.io/installation_guide.html). Once FELiCS is installed, the notebooks can be executed directly.

```{admonition} Disclaimer
These guides are not intended to be comprehensive manuals. Instead, they serve as a starting point for experimentation and learning.
They highlight core functionality and provide illustrative examples, but they do not cover every possible use case or advanced feature.
Users are encouraged to adapt the notebooks to their own projects and explore further.
```

Additional How-To guides are planned and will be added over time, expanding the coverage of FELiCS and offering deeper insights into its capabilities.
The current collection should be seen as a foundation upon which more detailed documentation will continue to grow.


### 1. First steps in scripting with FELiCS

In this section we will learn how to:
- import flow fields from h5 files into FELiCS objects and export them


```{toctree}
:maxdepth: 1

01_FirstSteps/01_HowToImportExportH5File.md
```

### 2. Derivatives and integral quantities

In this section we will learn how to:
- Calculate the L2-norm of a quantity over the whole domain
- Calculate the gradient of a scalar field in different coordinate systems, with symmetry and with spectral dimensions
- Calculate the vorticity of a 2D vector field


```{toctree}
:maxdepth: 1
02_CalculateL2Norm/02_CalculateL2Norm.md
02_CalculateGradient_cartesian/02_CalculateGradient_cartesian.md
02_CalculateGradient_cylindrical/02_CalculateGradient_cylindrical.md
02_CalculateVorticity/02_CalculateVorticity.md

```


### 3. Performing modal analysis in FELiCS

In this section we will learn how to:
- Perform a modal analysis on a given baseflow with a FEliCS config file

```{toctree}
:maxdepth: 1
03_PerformingModalAnalysisFELiCS/03_PerformingModalAnalysisFELiCS.md
```
### 4. Loading modes and structural sensitivity in FELiCS

In this section we will learn how to:
- Load previously computed resolvent or eigenmodes from disk using the FELiCS I/O module.

```{toctree}
:maxdepth: 1
04_LoadingModesAndStructuralSensitivity/04_LoadingModesAndStructuralSensitivity.md
```