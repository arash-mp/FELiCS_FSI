# MeshFilePath ([`Case`](index.md))

Path to the computational mesh file used for the FELiCS computation.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>None</code></dd>
</dl>

### Example
```json
"Case": {
    "MeshFilePath": "mesh.msh"
}
```

The mesh file defines the computational domain and discretization used by FELiCS.

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.

FELiCS currently supports Gmsh mesh files (`.msh`).