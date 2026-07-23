# m ([`Case`](index.md))

Defines the wavenumber associated with the spectral direction in FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>float</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>0.0</code></dd>
</dl>

### Example
```json
"Case": {
    "m": 1.0
}
```

This parameter is relevant for 2D base flows ([`nDim: 2`](nDim.md)) with one homogeneous direction represented spectrally in FELiCS.

Its interpretation depends on the selected [`CoordinateSystem`](CoordinateSystem.md):

- `"Cartesian"`:
  
  `m` represents the spanwise wavenumber associated with the homogeneous direction.

- `"Cylindrical"`:
  
  `m` represents the azimuthal mode number and should therefore take **integer** values.

This parameter is ignored for 3D base flows ([`nDim: 3`](nDim.md)).