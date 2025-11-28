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
autoapi_ignore = ["main.py", "runInputOutput.py", "runModal.py", "runResolvent.py"]

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

html_logo = "_static/logo_v2_full.png"

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
html_favicon = '_static/favicon.ico'