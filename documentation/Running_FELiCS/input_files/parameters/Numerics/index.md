# Numerics parameters

The `Numerics` category contains parameters related to the numerical discretization and solver configuration used by FELiCS.

These parameters control the polynomial approximation, eigenvalue targeting, and number of computed solutions.

## Contents

| Parameter | Description |
| - | - |
| [`EigenValueGuess`](EigenValueGuess.md) | Target eigenvalues used in modal analysis |
| [`nSolut`](nSolut.md) | Number of solutions computed by FELiCS |
| [`PolynomialOrder`](PolynomialOrder.md) | Polynomial order associated with each transported variable |

## Example
```json
"Numerics": {
    "nSolut": 10,
    "EigenValueGuess": ["1.0+0.1j"],
    "PolynomialOrder": {
        "u": 2,
        "T": 1
    }
}
```

## Related documentation

--

## Detailed content
```{toctree}
:maxdepth: 1

EigenValueGuess.md
nSolut.md
PolynomialOrder.md
```