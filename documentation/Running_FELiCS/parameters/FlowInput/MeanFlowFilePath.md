# MeanFlowFilePath ([`FlowInput`](index.md))

Path to the mean-flow input file used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"FlowInput": {
    "MeanFlowFilePath": "meanflow.fel"
}
```

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.

The expected file type depends on the value of [`needInterpolation`](../Case/needInterpolation.md):

| `needInterpolation` | Typical file type | Description |
| - | - | - |
| `true` | `.fel` | Mean flow defined on an external CFD grid and interpolated onto the FELiCS mesh |
| `false` | `MeanFlow.h5` | Mean flow previously exported by FELiCS on the same computational mesh |

See the [mean-flow input documentation](../../mean_files.md) for more details on supported formats.