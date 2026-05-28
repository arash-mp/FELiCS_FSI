# Reaction ([`Case`](index.md))

Enables reaction-related quantities in FELiCS.

<dl>
  <dt><strong>Type</strong></dt>
  <dd><em>boolean</em></dd>

  <dt><strong>Default</strong></dt>
  <dd><code>false</code></dd>
</dl>

### Example
```json
"Case": {
    "Reaction": true
}
```

This parameter is mainly relevant for reacting-flow simulations, telling the solver to load needed mean flow variables.