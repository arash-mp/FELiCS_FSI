.. FELiCS documentation master file, created by
   sphinx-quickstart on Wed Jun 21 16:35:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. .. automodule:: FELiCS
   :members:
   :no-undoc-members:

===============================
Welcome to FELiCS 
===============================

**FELiCS** (*Finite Element Linearized Combustion Solver*) is a Python-based CFD tool developed at the **Laboratory for Flow Instabilities and Dynamics**. It provides a powerful framework for analyzing and controlling fluid flows by solving the **linearized flow equations around a mean state**.

Originally built for academic research, FELiCS is designed to extend the application of linearized methods to real-world engineering problems. It supports a wide range of physical phenomena—including turbulence, heat and mass transport, chemical reactions, and acoustics—through a modular and extensible design.

FELiCS is:

- Built on top of **FEniCS** for finite element discretization
- Written in Python
- Integrated with **PETSc/SLEPc** for high-performance linear algebra
- Equipped for **linear stability analysis**, **resolvent analysis**, and **input-output analysis**
- Multi-physics: turbulence, mass transport, heat transport, chemical reactions, acoustics and many more
- supports low-Mach compressibitly
- Usable via the **command line** or integrated into **custom Python scripts**

This documentation will guide you through everything you need to know to get started with FELiCS: from installation to running analyses and diving into the math behind it.

--------------------------
Getting Started
--------------------------

1. 🚀 `Installation Guide <installation_guide.html>`_ 

   Step-by-step instructions for installing FELiCS and its dependencies.

2. ⚙️ `Running FELiCS <Running_FELiCS/index.html>`_

   Learn how to configure and run a simulation. This section also describes all available options and settings.

3. 🎓 `Tutorials <Tutorials/index.html>`_ 

   Experience modal and resolvent analysis.

--------------------------
Explore Further
--------------------------

* 📖 `Governing Equations <GoverningEquations/index.html>`_

  An in-depth look at the equations FELiCS solves and their physical meaning.

* 🛠️ `How-To-Guides <How-To-Guides/index.html>`_

  Additional tips and instructions for specific tasks and advanced features.

* 🧩 `API Documentation <autoapi/index.html>`_

  More information on classes, methods and functions for developers and users.

--------------------------
Content
--------------------------

.. toctree::
   :maxdepth: 1
   :glob:

   installation_guide
   Running_FELiCS/index
   Tutorials/index
   GoverningEquations/index
   How-To-Guides/index

   