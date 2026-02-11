# FELiCS Code Refactoring Guidelines

This document describes the **refactoring workflow for the FELiCS codebase**, with the goal of enforcing **PEP8 compliance**.

---

## Goals

- Enforce consistent **PEP8 naming and formatting**
- Improve readability and maintainability for new contributors
- Prevent cyclic imports and unnecessary heavy imports at startup

---

## 1. Prerequisites & Installation

Before beginning, configure your Python environment with the required tooling.

### Required Tools

```bash
pip install libcst isort black flake8 pylint rope
```

### Tool Overview

- **libcst**  
  Concrete Syntax Tree parser used to safely rewrite Python code while preserving comments and formatting.

- **isort**  
  Automatically sorts and groups imports (standard library, third-party, local).

- **black**  
  Enforces a consistent formatting style and line length (80 characters).

- **flake8 / pylint**  
  Static analysis tools for detecting style violations and potential logic errors.

- **rope**  
  Advanced refactoring library used for advanced rename operations.

---

## 2. Automated Refactoring Pipeline

We use a custom LibCST script, **`felics_style_rewriter.py`**, to perform safe, large-scale syntax refactoring.

The pipeline consists of two passes:

- **Pipeline A (`defs`)**: Refactors definitions and records renames
- **Pipeline B (`usages`)**: Applies recorded renames consistently across the codebase

---

## Step 1: Generate Definitions and Symbol Map

```bash
python code_refactoring/felics_style_rewriter.py defs ./src/
```

This step runs **Pipeline A (`defs`)** and performs the following actions on every Python file:

### What this step does

1. **Lifts imports to module scope**  
   - Moves nested `import` statements to the top of the file.
   - Improves readability and avoids hidden dependencies.

2. **Renames class definitions to PascalCase**  
   - Converts class names to PEP8-compliant PascalCase.
   - Preserves leading underscores.
   - Examples:
     - `meanFlowClass` → `MeanFlowClass`
     - `timeStepping` → `TimeStepping`

3. **Renames function definitions to snake_case**  
   - Converts function and method names to snake_case.
   - Examples:
     - `ComputeFlux` → `compute_flux`
     - `getPressureField` → `get_pressure_field`

4. **Normalizes trailing commas**  
   - Adds trailing commas to multi-argument calls and parameter lists.
   - Ensures compatibility with `black` and cleaner diffs.

5. **Generates `felics_renames.json`**  
   - Records all `(file, old_name, new_name)` triples.
   - This file is required for the second pass.

### What this step does NOT do

- Does not update usages or imports
- Does not rename files or directories
- Does not fix architectural or semantic issues

---

## Step 2: Review and Clean `felics_renames.json`

The `defs` step produces a file named **`felics_renames.json`**, which acts as the authoritative rename map.

### ⚠️ Manual Review Required

Some entries may be **false positives** and must be removed before proceeding.

### Common False Positives to Remove

Delete mappings such as:

- `t → T` (often represents time)
- `get_ksp → getKSP` (PETSc solvers)
- `Get_size → getSize` (external APIs)
- `Get_vecs` and similar backend-specific functions

Failing to clean this file may introduce incorrect renaming.

---

## Step 3: Apply Renames Across Usages

```bash
python code_refactoring/felics_style_rewriter.py usages ./src/
```

This runs **Pipeline B (`usages`)**, which:

- Updates all call sites to renamed functions and classes
- Fixes imports and aliases
- Renames attribute access consistently
- Reapplies trailing commas where needed

This step must only be run **after** reviewing `felics_renames.json`.

💡 Tip: If you are modifying or refactoring a single existing script, you can safely run the usages command on that file alone instead of the entire source tree:

```bash
python code_refactoring/felics_style_rewriter.py usages ./path/to/script.py
```
This is useful for incremental refactoring or when testing changes locally without touching unrelated files.
---

## 3. Manual Corrections

Automated refactoring handles most syntax changes, but some operations must be done manually.

---

## 3.1 File Renaming Strategy

All filenames must follow **snake_case**.

### Runner Scripts

Rename entry-point scripts:

- `run_input_output.py`
- `run_modal.py`
- `run_resolvent.py`

### Core Source Files

- `Parameters/Config.py`
- `Fields/FlunctuationClass.py`
- `Fields/FieldProperties.py`
- `Fields/MeanflowClass.py`

### Equations & Handlers

Standardize handler filenames:

- `Equations/dependentVariables/EnergyHandler.py`
- `Equation/dependentVariables/EquationOfStateHandler.py`
- `Equation/dependentVariables/HeatReleaseHandler.py`
- `Equation/dependentVariables/MomentumHandler.py`
- `Equation/dependentVariables/ReactionHandler.py`

Ensure all corresponding imports are updated.

---

## 3.2 Import Placement & Architecture Refinement

### 3.2.1 Lazy Loading Heavy Dependencies

Move heavy or optional imports inside the functions that use them.

#### `FELiCS/Misc/logging.py`
- Move `import IPython` inside `in_notebook()`
- Prevents crashes in environments without IPython

#### `Fields/MeanFlowClass.py`
- Move `from fenics import project` inside `calculate_species_enthalpy()`
- Avoids unnecessary FEniCS startup cost

---

### 3.2.2 Fixing Cyclic Imports

Break circular dependencies by deferring imports.

#### `Fields/FieldProperties.py`

- Move  
  ```python
  from FELiCS.Fields.MeanFlowClass import MeanFlowClass
  ```  
  inside `is_mean_flow()`

- Move  
  ```python
  from FELiCS.Fields.MeanFlowClass import MeanFlowVertexValues
  ```  
  inside `is_mean_flow_vertex_values_class()`

This prevents import-time failures.

---

## 3.3 Missing Dependencies After Cleanup

Automated cleanup may remove imports required dynamically.

### Runner Scripts

Files affected:
- `run_resolvent.py`
- `run_modal.py`
- `run_input_output.py`

Add explicitly:

```python
from FELiCS.SpaceDisc.FEMSpaces import FEMSpaces
```

inside the functions 

- `run_resolvent`
- `run_modal`
- `run_input_output`

---

## Final Notes

- Refactor incrementally and test frequently
- Commit automated and manual changes separately
- Favor explicit imports and clarity over convenience

Following this guide ensures the FELiCS codebase remains clean, consistent, and contributor-friendly.
