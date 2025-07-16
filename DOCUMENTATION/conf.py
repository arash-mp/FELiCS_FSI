# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# TODO: this should work to exclude the attributes from being shown, but it doesn't.
def skip_member(app, what, name, obj, skip, options):
    if what == 'attribute':
        return True
    return skip

def setup(app):
    app.connect("autodoc-skip-member", skip_member)

project = 'FELiCS2.0'
copyright = '2025, Flow Group'
author = 'Flow Group'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
# sys.path.insert(0, os.path.abspath('../Reactions/'))

extensions = ['sphinx.ext.coverage', 'sphinx.ext.napoleon', 'myst_parser', 'autoapi.extension','sphinxcontrib.mermaid']#, 'sphinx.ext.inheritance_diagram'] # 'sphinx.ext.autodoc', 
autoapi_dirs = ['../src/']
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autoapi_options = ['members', 'undoc-members', 'no-private-members', 'show-inheritance', 'show-module-summary', 'special-members', 'no-imported-members']  # 'private-members', 'members' TODO: check with whole group if special members or not
# autoapi_python_class_content = 'both' # renders __init__ docstring and class docstring together
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']


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

#Myst Parser settings
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_dmath_double_inline = True
