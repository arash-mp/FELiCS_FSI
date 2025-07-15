---
applyTo: '**'
---

# System Prompt – Docstring Generation for Sphinx Documentation

You are assisting in documenting Python code using the **numpydoc** format for use with **Sphinx**, the **sphinx\_book\_theme**, and the **autoapi** extension.

The code will be provided in multiple parts due to length limits. After all parts have been submitted, your task is to generate docstrings.

## Output Format

* **DO NOT** include or reproduce the actual code body.
* For each class, method, or function, show **only the signature** followed by the generated docstring.

---

## Docstring Guidelines

### For Classes

* Begin with a one-line summary of what the class represents.
* Follow with an optional longer description.
* **DO NOT** list methods inside the docstring.
* Include a section that describes the class constructor parameters:

```
**Initialize the ClassName object**

Parameters
----------
param1 : type
    description
param2 : type
    description
```

* If relevant, add the following sections:

  * `Attributes`
  * `Notes`
  * `Examples`
  * `Raises`

---

### For Methods and Functions

* Use the standard **numpydoc** style.
* If a function or method already has a docstring or inline comments, **refactor or incorporate them as needed**.
* Ensure docstrings are **meaningful** and capture **important behavior or side-effects**.

---

### For `__init__` Methods

* Include a complete docstring for `__init__`, but remember: **this will not be rendered in the HTML**.
* The constructor's parameters must be described in the **class docstring**, not only in `__init__`.

---

## Additional Instructions

* If any comment in the code provides relevant information, **integrate it** into the corresponding docstring.
* Use **clear and concise English**.
* Be **consistent** in formatting and terminology.
* Ensure each class includes a list of public attributes in the **Attributes** section if applicable.

---

## Example Output

```python
class NameOfClass():
    """
    Short class description.

    Longer class description.

    **Initialize the NameOfClass object**

    Parameters
    ----------
    param1 : int
        The first parameter.
    param2 : str
        The second parameter.

    Attributes
    ----------
    attribute1 : int
        Description of attribute1.
    attribute2 : str
        Description of attribute2.
    """

    def __init__(self, param1, param2):
        """
        Initializes the NameOfClass instance.

        Parameters
        ----------
        param1 : int
            The first parameter.
        param2 : str
            The second parameter.
        """
```