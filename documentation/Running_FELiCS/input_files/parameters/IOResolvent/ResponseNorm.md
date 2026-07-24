# ResponseNorm ([`IOResolvent`](index.md))

Defines the norm used to measure response amplitudes in resolvent analyses.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"TKE"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"TKE"` | Turbulent kinetic energy norm |
| `"Chu"` | Chu compressible energy norm |

### Example
```json
"IOResolvent": {
    "ResponseNorm": "Chu"
}
```

This parameter is only relevant for [resolvent analyses](../../../../GoverningEquations/analysisMode_resolvent.md).