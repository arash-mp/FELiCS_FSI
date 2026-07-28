# FlowInput parameters

The `FlowInput` category contains parameters related to the mean-flow input used by FELiCS.

These parameters control how the mean flow is imported and processed before the linear analysis.

## Contents

| Parameter | Description |
| - | - |
| [`MeanFlowFilePath`](MeanFlowFilePath.md) | Path to the mean-flow input file |

## Example
```json
"FlowInput": {
    "MeanFlowFilePath": "meanflow.fel"
}
```

## Related documentation

- See FELiCS input [mean-flow files](../../mean_files.md)

## Detailed content
```{toctree}
:maxdepth: 1

MeanFlowFilePath.md
```