# CalculateAdjoint ([`Case`](index.md))

Enables the computation of adjoint modes in modal analysis.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>boolean</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>false</code></dd>
</dl>

### Example
```json
"Case": {
    "CalculateAdjoint": true
}
```

This parameter is only relevant when [`"AnalysisMode"`](AnalysisMode.md) is set to `"Modal"`.

Note that FELiCS relies on a *discrete adjoint* framework. 