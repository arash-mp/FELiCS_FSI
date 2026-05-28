# nSolut ([`Numerics`](index.md))

Defines the number of solutions computed by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>integer</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>1</code></dd>
</dl>

### Example
```json
"Numerics": {
    "nSolut": 10
}
```

The interpretation of `nSolut` depends on the selected analysis type:

- `"Modal"`:
  
  Number of eigenmodes computed near each target eigenvalue defined by [`EigenValueGuess`](EigenValueGuess.md).

- `"Resolvent"` or `"InputOutput"`:
  
  Number of forcing/response mode pairs computed for each frequency defined by [`Omegas`](../IOResolvent/Omegas.md).

This parameter is passed internally to the PETSc/SLEPc eigensolvers and singular-value decomposition routines used by FELiCS.