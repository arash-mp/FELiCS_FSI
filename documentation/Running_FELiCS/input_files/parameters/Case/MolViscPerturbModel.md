# MolViscPerturbModel ([`Case`](index.md))

Defines the model used for fluctuating molecular viscosity in FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>dictionary</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>{"type":"Constant","Constants":{"Viscosity":1.0}}</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"Constant"` | Assume a frozen molecular viscosity |
| `"Sutherland mean"` | Computes fluctuating molecular viscosity using a linearized Sutherland-law model around the mean temperature |

### Example

The parameter should be provided as a dictionary containing:
- `"type"`: perturbation model name,
- `"Constants"`: dictionary containing model-specific constants.

```json
"Case": {
    "MolViscPerturbModel": {
        "type": "Sutherland mean",
        "Constants": {
            "Ts": 170.672
        }
    }
}
```

For `"Constant"`, the following constant is currently required:
- `"Viscosity"`: constant fluctuating viscosity value.

For `"Sutherland mean"`, the following constant is currently required:
- `"Ts"`: Sutherland temperature constant.