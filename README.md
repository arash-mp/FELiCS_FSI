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

## [Click here to access Documentation](https://felics2-0-laboratory-for-flow-instabilities-and--112f91add91c56.gitlab-pages.tu-berlin.de/)

## Guides
- [Installation Guide](DOCUMENTATION/installation_guide.md)
- [Tutorial 1](DOCUMENTATION/tutorial_1.md)
- [Governing Equations](https://www.overleaf.com/project/65fc3be353efc8edb72ec27e)
- [Best Practice Guidelines for Coding](https://www.overleaf.com/project/664f097a091050639ff51820)
- [A nice guide to Solving PDEs in Python](https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf)

## Validation Cases
The validation cases can be downloaded from [this gitlab repository] (https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics-tests>).

## How to write a guide or documentation file
- All guides and documentation files should be written in markdown or rst file format. 
- Here is a [template](DOCUMENTATION/markdown_template.md)
- To add a file to the documentation, put the markdown or rst file in the `docs/`folder
- Then, add the filename to the index file: `docs/index.rst`
- add a link to the file, here in the README file. This README will not be displayed on the web-based documentation. It is only ment for gitlab users, who are working on the code
- If you add content to this file (`felics2.0/README.md`) that shall also appear on the web-based version, you have to manually copy the content to 

## how to write docstrings
- We use the numpydoc docstring format
- see https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html for good examples
- see PostProcessing.py for example
- https://numpydoc.readthedocs.io/en/latest/format.html
- VSCode extension: autoDocstring - Python Docstring Generator is helpful!
