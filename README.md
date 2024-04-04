# FELiCS 2.0

```bash
(         (               (
)\ )      )\ )       (    )\ )
(()/(  (  (()/( (    )\   (()/(
/(_)) )\  /(_)))\  (((_)  /(_))
(_)_)((_) (_)) ((_) )\___ (_))
| __|| __|| |   (_)((/ __|/ __|
| _| | _| | |__ | | | (__ \__ \
|_|  |___||____||_|  \___||___/
```


## Guides and Documentation

The web-based documentation can be found here: [laboratory-for-flow-instabilities-and-dynamics.gitlab-pages.tu-berlin.de/felics2.0/](https://laboratory-for-flow-instabilities-and-dynamics.gitlab-pages.tu-berlin.de/felics2.0/)

## Guides
- [Installation Guide](DOCUMENTATION/source/installation_guide.md)
- [Tutorial 1](DOCUMENTATION/source/tutorial_1.md)
- [Best Practice Guidelines for Coding](https://www.overleaf.com/project/61dd7ece45e8fc038220b1b0)
- [A nice guide to Solving PDEs in Python](https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf)

## Validation Cases
The validation Cases can be found [here](https://tubcloud.tu-berlin.de/s/3MQCKgK7JKGSDdx)

## How to write a guide or documentation file
- All guides and documentation files should be written in markdown or rst file format. 
- Here is a [template](DOCUMENTATION/source/markdown_template.md)
- To add a file to the documentation, put the markdown or rst file in the `docs/source`folder
- Then, add the filename to the index file: `docs/source/index.rst`
- add a link to the file, here in the README file. This README will not be displayed on the web-based documentation. It is only ment for gitlab users, who are working on the code
- If you add content to this file (`felics2.0/README.md`) that shall also appear on the web-based version, you have to manually copy the content to 

## how to write docstrings
- We use the numpydoc docstring format
- see https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html for good examples
- see PostProcessing.py for example
- https://numpydoc.readthedocs.io/en/latest/format.html
- VSCode extension: autoDocstring - Python Docstring Generator is helpful!