# How-To-Guides

Welcome to How-To-Guides!

**TODO: put the following bullet points into text:**
- these markdowns are based on jupyter notebooks
- all the jupyter notebooks and the files that they need are in the folder "how_tos" (link here!) 
- remember to install FELiCS as a package before use (link the intallation guide here)
- these how-to guides are by no means complete but a starting point for exploring how to use FELiCS as a python package
- further How-To Guides are planned


### First steps in scripting with FELiCS

In this section we will learn how to:
- import flow fields from h5 files into FELiCS objects and export them


```{toctree}
:maxdepth: 1

01_ImportExport.md
```

### Derivatives and integral quantities

In this section we will learn how to:
- Calculate the L2-norm of a quantity over the whole domain
- Calculate the gradient of a scalar field in different coordinate systems, with symmetry and with spectral dimensions
- Calculate the vorticity of a 2D vector field


```{toctree}
:maxdepth: 1
02_CalculateL2Norm.md
02_CalculateGradient_cartesian.md
02_CalculateGradient_cylindrical.md
02_CalculateVorticity.md

```


### Performing modal analysis in FELiCS

In this section we will learn how to:
- How to perform a modal analysis on a given baseflow with a FEliCS config file

```{toctree}
:maxdepth: 1
03_GeneralEigenValueProblem.md
```
