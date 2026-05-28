# IOResolvent parameters

The `IOResolvent` category contains parameters related to resolvent and input/output (IO) analyses performed by FELiCS.

These parameters define the forcing configuration, response norms, and analyzed frequencies.

## Contents

| Parameter | Description |
| - | - |
| [`ForcingBoundaryIndices`](ForcingBoundaryIndices.md) | Boundary indices where forcing is applied in IO analysis |
| [`ForcingMode`](ForcingMode.md) | Type of forcing applied in the analysis |
| [`ForcingNorm`](ForcingNorm.md) | Norm used for the forcing term |
| [`ResponseNorm`](ResponseNorm.md) | Norm used for the response term |
| [`Omegas`](Omegas.md) | Angular frequencies analyzed of the run |

## Example
```json
"IOResolvent": {
    "Omegas": [1.0, 2.0],
    "ForcingMode": "Body",
    "ForcingNorm": "Chu",
    "ResponseNorm": "Chu"
}
```

## Related documentation

See the documentation pages about the [resolvent](../../../GoverningEquations/analysisMode_resolvent.md) and [input/output](../../../GoverningEquations/analysisMode_input_output.md) analyses.

## Detailed content
```{toctree}
:maxdepth: 1

ForcingBoundaryIndices.md
ForcingMode.md
ForcingNorm.md
ResponseNorm.md
Omegas.md
```