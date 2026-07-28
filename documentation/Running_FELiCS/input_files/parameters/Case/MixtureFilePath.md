# MixtureFilePath ([`Case`](index.md))

Path to the mixture-definition file used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"Case": {
    "MixtureFilePath": "mixture.json"
}
```

The mixture file defines the thermodynamic, transport, and chemical properties of the fluid mixture considered in the simulation.

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.

See the [mixture-file documentation](../../mix_files.md) for more details on the expected file format.