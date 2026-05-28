# Omegas ([`IOResolvent`](index.md))

Defines the angular frequencies analyzed by FELiCS in resolvent or input/output analyses.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>list of floats or complex values</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>[]</code></dd>
</dl>

### Example
```json
"IOResolvent": {
    "Omegas": [1.0, 2.0, 5.0]
}
```

Complex frequencies should be provided as strings:

```json
"IOResolvent": {
    "Omegas": ["1.0+0.1j", "2.0-0.2j"]
}
```

Real and complex frequencies may also be combined:

```json
"IOResolvent": {
    "Omegas": [1.0, "2.0+0.1j"]
}
```

Each entry defines one angular frequency analyzed by FELiCS.

This parameter is not used when [AnalysisMode](../Case/AnalysisMode.md) is set to `Modal`.