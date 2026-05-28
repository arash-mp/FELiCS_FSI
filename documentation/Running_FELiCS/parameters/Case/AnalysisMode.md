# AnalysisMode ([`Case`](index.md))

Controls the type of linear analysis performed by the FELiCS run.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"Modal"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"Modal"` | Indicates a linear modal (stability) analysis |
| `"Resolvent"` | Indicates a resolvent (optimal input/output) analysis |
| `"InputOutput"` | Indicates an input/output (linear response or harmonic response) analysis |

### Example
```json 
"Case": { 
    "AnalysisMode": "Resolvent"
}
```

For more information about the theory behind each type of analysis, see the pages about [modal](../../../GoverningEquations/analysisMode_modal.md), [input-output](../../../GoverningEquations/analysisMode_input_output.md) and [resolvent](../../../GoverningEquations/analysisMode_resolvent.md) analyses.