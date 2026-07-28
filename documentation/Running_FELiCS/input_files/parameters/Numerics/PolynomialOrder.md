# PolynomialOrder ([`Numerics`](index.md))

Defines the polynomial interpolation order associated with each transported variable in FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>dictionary</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"Numerics": {
    "PolynomialOrder": {
        "u": 2,
        "T": 1,
        "rho": 1
    }
}
```

The parameter should be provided as a dictionary whose entries associate:
- a transported variable,
- with its corresponding polynomial interpolation order.

The variable names should be consistent with those defined in [`SetOfEquations`](../Case/SetOfEquations.md).