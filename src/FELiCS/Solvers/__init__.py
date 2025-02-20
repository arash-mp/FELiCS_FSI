"""
Resolvent Operator Package
==========================

This package provides tools for working with resolvent operators in the context of singular value decomposition (SVD) and linear equation systems, leveraging the PETSc library for efficient matrix and vector operations.

Main Features
-------------
- **ResolventOperator**:
  A matrix-free representation of the resolvent operator multiplied by its Hermitian transpose, supporting singular value decomposition.

- **Linear Equation Solvers**:
  Utilities for creating and solving linear systems efficiently using preconditioned solvers with PETSc.

Modules
-------
- `resolvent_operator`:
  Contains the `ResolventOperator` class and helper functions for creating and solving linear systems.

Dependencies
------------
This package relies on the `petsc4py` library, which provides Python bindings for PETSc.

Usage
-----
Import the main components as follows:

    from resolvent_operator import ResolventOperator, createEquationSystemSolver, solveEquationSystemWithPredefinedSolver

Example:

    # Create a solver
    solver = createEquationSystemSolver(A)

    # Solve a linear system
    solution = solveEquationSystemWithPredefinedSolver(solver, b)

    # Work with the ResolventOperator
    operator = ResolventOperator(ResolventMatrix, W_FEM, W_forcing, W_response, Restrictor_forcing, Restrictor_response)
    size = operator.getSize()

"""
