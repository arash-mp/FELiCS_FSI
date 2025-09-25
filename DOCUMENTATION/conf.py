project = 'FELiCS2.0'
copyright = '2025, Flow Group'
author = 'Flow Group'

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

extensions = ['sphinx.ext.coverage', 'sphinx.ext.napoleon', 'myst_parser', 'autoapi.extension','sphinxcontrib.mermaid']
autoapi_dirs = ['../src/']
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autoapi_options = ['members', 'undoc-members', 'no-private-members', 'show-inheritance', 'show-module-summary', 'special-members', 'no-imported-members']

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

html_logo = "_static/logo.png"

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Myst Parser settings
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_dmath_double_inline = True


# Where the switcher JSON lives (root of each deployment)
html_js_files = ['version-switcher.js']
_SWITCHER_JSON = "/versions.json"

# Figure out which version we’re currently viewing:
# In your Pages parallel deploys, tags render under /<tag>/ and “latest” at /
_version = os.environ.get("CI_COMMIT_TAG", "latest")

html_theme_options = {
    "navbar_end": ["version-switcher"],
    "switcher": {
        "json_url": "/versions.json",
        "version_match": os.environ.get("CI_COMMIT_TAG", "latest"),
    },
}
