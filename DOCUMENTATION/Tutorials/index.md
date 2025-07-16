# Tutorials

Welcome to the FELiCS Tutorials!

These step-by-step guides are designed to help you get started with FELiCS, from basic case setup and mesh generation to running modal and resolvent analyses. Each tutorial builds on the previous one, introducing new concepts and tools.

Below you will find an overview of the available tutorials, with links and a short description for each.

```{admonition} Prerequisites
Before starting these tutorials, ensure that FELiCS is installed on your system.  
If you haven't done so yet, see the [Installation Guide](../installation_guide.md) for step-by-step instructions.
```

## [Tutorial 1: Cylinder Wake – Solve the Base Flow](cylinder_wake.md)
Learn how to set up your first FELiCS case by computing the 2D base flow around a cylinder at Reynolds number 50.

You will:
- Define a case in FELiCS.
- Create a case folder with a mesh and base flow files.
- Run the base flow solver and visualize results in Paraview.

---

## [Tutorial 2: Modal Analysis](modal_analysis.md)
Perform a linear stability (eigenvalue) analysis of the base flow computed in Tutorial 1.

You will:
- Define and apply boundary conditions.
- Set up and run a modal analysis case.
- Postprocess and visualize eigenmodes and spectra.


---

## [Tutorial 3: Resolvent Analysis](Resolvent_Analysis.md)
Explore frequency response and amplification mechanisms in a constricted pipe (stenosis) using resolvent analysis.

You will:
- Generate a mesh with GMSH.
- Import and process mean flow data.
- Set up and run resolvent analysis for different frequencies.
- Postprocess and plot gains and mode shapes.

---

## [Tutorial 4: Input-Output Analysis](input_output_analysis.md)
Lorem Ipsum
You will:
- Do stuff

---

## Additional Resources
- [GMSH Quick Reference](gmsh.md):  Handy guide to the most important GMSH Python API commands for mesh generation.

---

Each tutorial folder contains all necessary scripts, input files, and example outputs.  
For further details on settings and file formats, see the [FELiCS settings documentation](../Running_FELiCS/FELiCS_settings.md).


```{toctree}
:maxdepth: 1

cylinder_wake.md
modal_analysis.md
Resolvent_Analysis.md
input_output_analysis.md
```
