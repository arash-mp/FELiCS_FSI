# Tips and tricks for creating documentations like this


## Required Packages
1) sphinx: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/install.html 
1b) probably run `sudo apt-get install python3-sphinx --fix-missing`
2) Book Theme: https://sphinx-themes.org/sample-sites/sphinx-book-theme/ 
3) myst-parser: https://myst-parser.readthedocs.io/en/v0.17.1/sphinx/intro.html, `pip install myst-parser`

# commands

- create dummy folder: navigate to `project-name/docs` folder, run `sphinxs-quickstart`
- generate html by running `make html` in `docs` folder
-  

# API generation (also works but is more manual)
`sphinx-apidoc -o ./source/src ../src`

# AutoAPI (works for entire project folder)
https://sphinx-autoapi.readthedocs.io/en/latest/tutorials.html

`sphinx-build -b html . _build`




