# Solution files

FELiCS exports computed mode shapes and harmonic responses in individual `.h5` files.

These files are produced by:
- `"Modal"` analyses,
- `"Resolvent"` analyses,
- `"InputOutput"` analyses.

chosen with [`Case.AnalysisMode`](../input_files/parameters/Case/AnalysisMode.md).

## Contents

Solution files contain complex fields associated with the computed solution:

| Analysis type | File contents |
| - | - |
| `"Modal"` | Eigenmodes associated with the computed eigenvalues |
| `"Resolvent"` | Forcing and response singular vectors associated with the computed gains |
| `"InputOutput"` | Complex harmonic response fields |

The exported fields correspond to the fluctuating variables of the FELiCS **state vector**, as defined by [`SetOfEquations`](../input_files/parameters/Case/SetOfEquations.md).

The fields are exported on the FELiCS [output mesh](FELiCSmesh_files.md).

Complex fields are stored by separating their real (`variable_real`) and imaginary (`variable_imag`) parts.

## File naming

FELiCS generates solution file names from:
- the analysis type,
- the mode type,
- the angular frequency or eigenvalue,
- and, for resolvent analyses, the gain number.

Typical examples are:

| Analysis type | Example file name |
| - | - |
| `"Modal"` | `Mode_Modal_Adjoint_Omega_0.727-0.001j.h5` |
| `"Resolvent"` | `Mode_Resolvent_Forcing_Omega_0.700_GainNb_0.h5` |
| `"InputOutput"` | `Mode_Input_output_Response_Omega_628.319.h5` |

## Attributes

Solution files include metadata attributes describing the exported mode.

Typical attributes include:

| Attribute | Description |
| - | - |
| `omega` | Angular frequency or eigenvalue associated with the mode |
| `m` | Wavenumber associated with the spectral direction |
| `gain` | Gain associated with the mode, for resolvent and input/output analyses |
| `number` | Gain number, for resolvent analyses |