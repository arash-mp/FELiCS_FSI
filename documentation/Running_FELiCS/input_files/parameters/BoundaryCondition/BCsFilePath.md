# BCsFilePath ([`BoundaryCondition`](index.md))

Path to the boundary-condition file used by FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"BoundaryCondition": {
    "BCsFilePath": "boundaries.json"
}
```

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.