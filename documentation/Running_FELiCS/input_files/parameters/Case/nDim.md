# nDim ([`Case`](index.md))

Defines the number of spatial dimensions resolved by the computational mesh.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>integer</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>2</code></dd>
</dl>

| Accepted values | Description |
| - | - |
| `2` | Two-dimensional computational mesh |
| `3` | Three-dimensional computational mesh |

### Example
```json
"Case": {
    "nDim": 2
}
```

For `2D` simulations, FELiCS may optionally represent an additional homogeneous direction spectrally through the [`m`](m.md) parameter.

For `3D` simulations, the [`m`](m.md) parameter is ignored.