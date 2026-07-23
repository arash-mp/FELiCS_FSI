# Output files

After running an analysis, FELiCS generates several output files containing:
- computed solutions (`.h5`),
- interpolated mean-flow fields (`.h5`),
- computational meshes (`.h5`),
- solution spectra (`.csv`),
- and execution logs (`.log`).

For visualization purposes, FELiCS also automatically generates associated `.xmf` metadata files compatible with ParaView.

The generated files depend on the selected analysis type and export options.

The different output files are described in the following sections:

| File | Description |
| - | - |
| [Solution files](solution_files.md) | Modal, resolvent, or input/output solutions exported by FELiCS |
| [FELiCS Mean-flow files](FELiCSmeanflow_files.md) | Mean-flow fields exported onto the FELiCS mesh |
| [FELiCS Mesh files](FELiCSmesh_files.md) | Exported FEM meshes used by FELiCS |
| [Spectrum files](spectrum_files.md) | Eigenvalue or gain spectra exported as CSV files |
| [Log files](log_files.md) | Solver logs and execution information |

```{toctree}
:maxdepth: 1

solution_files
FELiCSmeanflow_files
FELiCSmesh_files
spectrum_files
log_files
```