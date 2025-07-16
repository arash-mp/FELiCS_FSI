# How-To-Guides

Welcome to How-To-Guides!

**TODO: re-format this index site and add some text**


### First steps in scripting with FELiCS

In this section we will learn how to:
- import flow fields from h5 files into FELiCS objects and export them


```{toctree}
:maxdepth: 1

01_ImportExport.md
```

### Derivatives and integral quantities

In this section we will learn how to:
- Calculate the gradient of a scalar field in different coordinate systems, with symmetry and with spectral dimensions
- Calculate the vorticity of a 2D vector field
- Calculate the L2-norm of a quantity over the whole domain

```{toctree}
:maxdepth: 1
02_CalculateGradient_cartesian.md
02_CalculateGradient_cylindrical.md
02_CalculateVorticity.md
02_CalculateL2Norm.md
```


### Performing modal analysis in FELiCS

In this section we will learn how to:
- How to perform a modal analysis on a given baseflow
- How to calculate the sensitivity of a eigenvalue with respect to a given baseflow

```{toctree}
:maxdepth: 1
03_GeneralEigenValueProblem.md
```

### In progress 
- Sensitivity tutorials: how to calculate eigenvalue sensitivity w.r.t. baseflow; how to calculate structural sensitivity
- tensorial framework tutorial


### Ideas for the future
- How to add custom equations
- how to import SPOD mode into felicsfunction and compute alignment with resolvent mode
- calculate pseudospectrum of a linear operator
