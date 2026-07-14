# ForcingMode ([`IOResolvent`](index.md))

Defines how forcing is applied in resolvent or input/output analyses.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"Body"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"Body"` | Applies forcing within the computational domain |
| `"Boundary"` | Applies forcing on selected boundaries |

### Example
```json
"IOResolvent": {
    "ForcingMode": "Boundary",
    "ForcingBoundaryIndices": [1, 3]
}
```

When `"Boundary"` forcing is selected, the forced boundaries should be specified through [`ForcingBoundaryIndices`](ForcingBoundaryIndices.md).

This parameter is only relevant when [Case.AnalysisMode](../Case/AnalysisMode.md) is set to `"InputOutput"`.