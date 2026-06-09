# SetOfEquations ([`Case`](index.md))

Defines the set of equations and transported variables used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>dictionary</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

**Available equations:**
| Equation Type | Form | Description |
| - | - | - |
| `"Mass"` | `"Continuity"` | [Mass conservation](../../../../GoverningEquations/equation_momentum_mass_species_energy.md#mass-equation) |
| `"Momentum"` | `"NSPrimitive"` | [Navier–Stokes momentum equations](../../../../GoverningEquations/equation_navier_stokes.md) |
| `"Species"` | `"Non-conservative"` | [Species transport equation](../../../../GoverningEquations/equation_species.md) |
| `"Energy"` | `"Enthalpy"` | Enthalpy transport equation (**DEPRECATED**) |
| `"Energy"` | `"primitive-p"` | [Pressure transport equation](../../../../GoverningEquations/equation_energy.md#1-energy-pressure-equation) |
| `"Energy"` | `"ProgressVariableLinear"` | [Scalar transport equation](../../../../GoverningEquations/equation_energy.md#2-progress-variable-equation) |
| `"EquationOfState"` | `"IdealGas"` | [Ideal gas equation of state](../../../../GoverningEquations/equation_state.md#1-ideal-gas-equation) |
| `"EquationOfState"` | `"Low-Mach"` | [Low-Mach approximation of equation of state](../../../../GoverningEquations/equation_state.md#2-low-mach-equation) |

### Example
```json
"Case": {
    "SetOfEquations": {
        "Momentum": {
            "Equation": "NSPrimitive",
            "Variable": "u"
        },
        "Energy": {
            "Equation": "primitive-p",
            "Variable": "T"
        }
    }
}
```

The parameter is provided as a dictionary whose entries define:
- the equation type,
- the specific form (in cases where several variations exist),
- the associated transported variable.

Each equation sub-block contains:
- `"Equation"`: equation model used by FELiCS,
- `"Variable"`: corresponding transported variable name.

### Notes

The variables defined in `SetOfEquations` should be consistent with other configuration parameters such as:
- [`PolynomialOrder`](../Numerics/PolynomialOrder.md)