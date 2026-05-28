# EigenValueGuess ([`Numerics`](index.md))

Defines the target eigenvalues used by FELiCS in modal analyses.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>list of floats or complex values</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>[]</code></dd>
</dl>

### Example
```json
"Numerics": {
    "EigenValueGuess": ["0.1+1.2j", "0.1-1.2j"]
}
```

Real eigenvalue targets may also be provided:

```json
"Numerics": {
    "EigenValueGuess": [0.0, 1.0]
}
```

Complex eigenvalues should be provided as strings.

Real and complex targets may also be combined:

```json
"Numerics": {
    "EigenValueGuess": [0.0, "0.1+1.2j"]
}
```

Each entry defines one eigenvalue target used internally by the eigensolver. This parameter is then passed to the `setTarget()` method of the SLEPc eigenvalue solver.

This parameter is only relevant when [Case.AnalysisMode](../Case/AnalysisMode.md) is set to `Modal`.