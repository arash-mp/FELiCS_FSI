# Export parameters

The `Export` category contains parameters controlling the export of FELiCS outputs and post-processing files.

## Contents

| Parameter | Description |
| - | - |
| [`ExportFolder`](ExportFolder.md) | Directory where FELiCS outputs are exported |

## Example
```json
"Export": {
    "ExportFolder": "results/",
}
```

## Related documentation

See the documentation on the FELiCS [output files](../../../index.md#output-files) to know which outputs are produced by the runs.

## Detailed content
```{toctree}
:maxdepth: 1

ExportFolder.md
```