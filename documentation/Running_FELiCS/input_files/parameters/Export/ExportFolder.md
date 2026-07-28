# ExportFolder ([`Export`](index.md))

Defines the directory where FELiCS exports output files.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>string</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>"./"</code></dd>
</dl>

### Example
```json
"Export": {
    "ExportFolder": "results/"
}
```

The path can be given either as:
- a relative path from the FELiCS running directory,
- or an absolute path.

FELiCS writes simulation outputs and post-processing files to this directory.