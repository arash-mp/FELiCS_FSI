# CoordinateSystem ([`Case`](index.md))

Defines the coordinate system used for the FELiCS computation.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"Cartesian"</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `"Cartesian"` | Cartesian coordinate system |
| `"Cylindrical"` | Cylindrical coordinate system |

### Example
```json
"Case": {
    "CoordinateSystem": "Cartesian"
}
```

The selected coordinate system affects the form of the linearized equations assembled by FELiCS.