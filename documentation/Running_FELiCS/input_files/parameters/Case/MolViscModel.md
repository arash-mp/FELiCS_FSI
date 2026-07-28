# MolViscModel ([`Case`](index.md))

Defines the **mean flow** molecular viscosity model used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"Constant"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"Constant"` | Uses a spatially uniform molecular viscosity defined by [`MolVisc`](MolVisc.md) |
| `"File"` | Imports the molecular viscosity field from the input [mean-flow file](../../mean_files.md) |

### Example
```json
"Case": {
    "MolViscModel": "Constant",
    "MolVisc": 1.5e-5
}
```

The selected viscosity model determines how the **mean flow** molecular viscosity field is constructed internally by FELiCS.