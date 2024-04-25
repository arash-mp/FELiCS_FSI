# Documentation Creation Guide for Developers

This guide explains how you can build the documentation locally to see how your docstrings and guides will render on the automatically created documentation html. 

## Installation
- `conda create --name felics_documentation_env`
- `conda activate felics_documentation_env`
- `conda install -c conda-forge sphinx`
- `conda install -c conda-forge myst-parser`
- `conda install -c conda-forge sphinx-book-theme`
- `conda install -c conda-forge sphinx-autoapi`

## locally build documentaion
- navigate to DOCUMENTATION folder `cd DOCUMENTATION`
- `sphinx-build -b html . _build`


## other commands
- create dummy folder: navigate to `project-name/docs` folder, run `sphinxs-quickstart`
- generate html by running `make html` in `docs` folder
- API generation (also works but is more manual) `sphinx-apidoc -o ./source/src ../src`



