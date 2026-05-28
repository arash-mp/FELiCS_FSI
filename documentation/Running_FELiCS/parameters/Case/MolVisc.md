# MolVisc ([`Case`](index.md))

Defines the molecular viscosity used by FELiCS when a constant viscosity model is selected.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>float</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>0.0</code></dd>
</dl>

### Example
```json
"Case": {
    "MolViscModel": "Constant",
    "MolVisc": 1.5e-5
}
```

The interpretation of `MolVisc` depends on the governing equations used in the simulation:

- incompressible equations:
  
  `MolVisc` represents the kinematic viscosity.

- compressible equations:
  
  `MolVisc` represents the dynamic viscosity.

This parameter is only used when [`MolViscModel`](MolViscModel.md) is set to `"Constant"`.