# ForcingBoundaryIndices ([`IOResolvent`](index.md))

Defines the indices of the boundaries where forcing is applied in input/output analyses.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>list of integers</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>[]</code></dd>
</dl>

### Example
```json
"IOResolvent": {
    "ForcingBoundaryIndices": [1, 3]
}
```

The indices should correspond to the physical boundary tags defined in the [computational mesh](../../gmsh_files.md).

This parameter is only used for `"InputOutput"` analyses with [`ForcingMode`](ForcingMode.md) set to `"Boundary"`. It is ignored for resolvent and modal analyses.

See also:
- [`ForcingMode`](ForcingMode.md)
- [`AnalysisMode`](../Case/AnalysisMode.md)