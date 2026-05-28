# SpeciesFilePath ([`Case`](index.md))

Path to the species-definition file used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"Case": {
    "SpeciesFilePath": "species.json"
}
```

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.

> 🚧 This parameter is currently not used by FELiCS, but is retained for compatibility with upcoming developments of the mixture initialization, relevant for reacting-flow cases.