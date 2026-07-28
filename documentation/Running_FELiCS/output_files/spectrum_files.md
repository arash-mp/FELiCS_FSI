# Spectrum files

FELiCS exports spectral solution information using `.csv` files.

These files are produced only for:
- `"Modal"` analyses,
- `"Resolvent"` analyses,

chosen with [`Case.AnalysisMode`](../input_files/parameters/Case/AnalysisMode.md).

## Modal spectrum

Modal analyses export the eigenvalue spectrum to:

```text
spectrum.csv
```

Typical columns include:

| Column | Description |
| - | - |
| `omega_direct_r` | Real part of the direct eigenvalue |
| `omega_direct_i` | Imaginary part of the direct eigenvalue |
| `omega_adjoint_r` | Real part of the adjoint eigenvalue |
| `omega_adjoint_i` | Imaginary part of the adjoint eigenvalue |

Adjoint columns are only present when [`Case.CalculateAdjoint`](../input_files/parameters/Case/CalculateAdjoint.md) is enabled.

### Example

```text
omega_direct_r,omega_direct_i,omega_adjoint_r,omega_adjoint_i
0.727343,0.001271,0.727343,-0.001271
```

## Resolvent spectrum

Resolvent analyses export the computed gains to:

```text
gains.csv
```

Typical columns include:

| Column | Description |
| - | - |
| `omega` | Angular frequency |
| `gain_0` | Leading resolvent gain |
| `gain_1` | Second resolvent gain |
| ... | Additional subleading gains |

### Example

```text
omega,gain_0,gain_1
0.7,1.100880e+07,5057.819664
0.8,1.394340e+07,2557.409904
```