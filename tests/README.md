# Tests for the FELiCS Package

This directory contains unit tests for individual methods in the **FELiCS** package.
The tests make use of very small grids and FEM functions that are initialized with random variables.
The folder structure inside the `tests` directory mirrors the (sub-)packages of FELiCS.

The main purpose of these tests is to cover functionality that is **not sufficiently tested** in the external repository **felics-tests**.
This applies especially to features required when **writing scripts with FELiCS**.

## Important Notes for Contributors

All contributors to FELiCS must ensure that **every new feature can be tested appropriately**.

You have two options:

1. Create or update a test in the external [felics-tests](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics-tests) repository so that your new feature is fully covered, **or**
2. Add sufficient **unit tests in this directory**.

No new feature should be merged without proper test coverage.

## Test Framework

We use the Python package **pytest**.

### Running Tests Locally

To run the tests on your local machine:

1. Install the `pytest` package (e.g., via `conda install pytest` or `pip install pytest`).
2. Install FELiCS as a Python package (see the installation guide).
3. Navigate into the `tests` directory.
4. Execute:

   ```bash
   python3 -m pytest
   ```

All test methods will be automatically discovered and executed.

### Automated Testing via GitLab

All tests in this folder are automatically executed by the **GitLab runner** whenever changes are pushed to the `development` branch.
You can view the success status of these tests via the **badge** in the main FELiCS README.

Before pushing to the `development` branch, always ensure that **all unit tests pass locally**.
Refer to the Wiki for the [Guidelines on merging into development](https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/wikis/Workflow-merging-into-development).

