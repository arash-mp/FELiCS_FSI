# ForcingNorm ([`IOResolvent`](index.md))

Defines the norm used to measure forcing amplitudes in resolvent analyses.

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
    "ForcingNorm": "Chu"
}
```

This parameter is only relevant for [resolvent analyses](../../../../GoverningEquations/analysisMode_resolvent.md).