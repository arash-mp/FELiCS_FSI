# TurbulenceModel ([`Case`](index.md))

Defines the mean flow turbulence model used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"None"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"None"` | No turbulence model |
| `"Constant"` | Uses a constant turbulent viscosity |
| `"File"` | Imports the turbulent viscosity field from the mean-flow input file |

### Example
```json
"Case": {
    "TurbulenceModel": "File"
}
```

The selected turbulence model determines how the mean flow turbulent viscosity field is constructed internally by FELiCS.