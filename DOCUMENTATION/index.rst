.. FELiCS2.0 documentation master file, created by
   sphinx-quickstart on Wed Jun 21 16:35:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. .. automodule:: Equations.Momentum.addMomentumEq
..     :members:
=====================================
FELiCS2.0
=====================================

The goal of the FELiCS project is to (further) develop a code (FELiCS) that applies linear analysis to multi-physics flow problems.
The software should provide a simple access to the related methods (Stability analysis, Resolvent analysis, Input-Output analysis) and allow a user to apply the code without significant knowledge about its details and implementation to flows in complex geometries.


`Guides for developers can be found in the Wiki <https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/wikis/home>`_

---------------
Content
---------------
.. toctree::
   :maxdepth: 1
   :glob:

   installation_guide
   tutorial_1
   documentation_creation


Validation Cases
-------------
The validation cases can be downloaded from `this gitlab repository <https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics-tests>`_ . 

Governing Equations
-------------
The Governing equations are documented in the following `overleaf document <https://www.overleaf.com/project/65fc3be353efc8edb72ec27e>`_

Coding Guidelines
-------------
Our Best Practice Guidelines for Coding in FELiCS are documented in the following `overleaf document <https://www.overleaf.com/project/664f097a091050639ff51820>`_
A nice guide to Solving PDEs in Python can be found on  `fenicsproject.org <https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf>`_


How to write a guide or documentation file
-------------

- All guides and documentation files should be written in markdown or rst file format. 
- The template is in docs/source/markdown_template.md
- To add a file to the documentation, put the markdown or rst file in the docs/source folder
- Then, add the filename to the index file: docs/source/index.rst

How to write docstrings
-------------
- We use the numpydoc docstring format
- see https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html for good examples
- see PostProcessing.py for example
- https://numpydoc.readthedocs.io/en/latest/format.html
- VSCode extension: autoDocstring - Python Docstring Generator is helpful!
=======
   Explanations/index
   Governing Equations/index
   Tutorials/index
   How-To-Guides/index

Indices and tables
-------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
