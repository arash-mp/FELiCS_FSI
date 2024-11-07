# Best coding practices
## Code Layout
Below are some general rules for FELiCS developers:
1. Indentation is done with four spaces and four spaces only.
2. For all variables, functions and classes the camelCase notation is used.
3. Only the English language is used in descriptive names and comments.
4. Long, descriptive names for classes, functions and variables (e.g. ''numberGridPoints'') are preferred to short cryptic ones (e.g. ''n'').
5. The maximum line length should not exceed 80 characters. If longer lines cannot be avoided, continuation lines should be used (see below).
6. There is a space character before and after every operator (''='', ''+'', ''-'', ''*'', ''/'', ''and'', ''or'', ''is'', ''$>$'', ''$<$'', ...)
7. There is no space before or after brackets (''( )'',''[ ]'' and ''\{ \}''), except they are preceded or followed by an operator.
8. There is no space before every comma and column, while there is one space after every comma and colon, except it is followed by a bracket.
9. Use ''redundant'' brackets and trailing commas to improve git version tracking.

## Continuation lines
When functions are defined or called, continuation lines should be used. The line to be continued ends with the open bracket. The continuing lines are preceded by a double hanging indent and contain one passed variable only. The closing bracket is on an own line with a double hanging indent.

**WRONG**: 
```python
def long_function_name1(var_one, var_two, 
    var_three, var_four):
    print(var_one)
foo = long_function_name1(var_one, var_two, var_three, var_four)
```
**CORRECT:**
```python
def long_function_name1(
        var_one, 
        var_two, 
        var_three, 
        var_four,
        ):
    print(var_one)
    
foo = long_function_name1(
    var_one, 
    var_two, 
    var_three, 
    var_four,
    )
```
## Imports
The following are rules related to imports in FELiCS:
1. Imports are performed at the beginning of every module/file. Use a new line for every new library.
2. Only import an entire library if you need more than ten classes or functions. If you need ten or less classes or functions use an absolute import.
3. If absolute import is used, begin a new file for every imported method or class, beginning with the second. Use a hanging intendation for the second and following methods and classes to highlight the absolute import. Furthermore, use the trailing commas to enhance git version tracking.
4. Both the libraries and the methods and classes are to be ordered alphabetically.
5. Absolute imports via * are discouraged as there is the risk of superseding other functions, variables and methods. The only exception is for the fenics library if more than ten function and classes are used.
6. The import section at the beginning of every module is divided in three parts: First the standard libraries, then the third party libraries, at last the local modules or libraries. These three sections are to be divided by a comments as shown in the example.

**WRONG**: 
```python
from local_library1 import method1, method2, class1
from fenics import * # Only the Function and the TrialFunction classes are needed
import sys, os
from local_library2 import *
```
**CORRECT:**
```python
# Standard libraries
import os
import sys

# Third party libraries
from fenics import (
    Function, 
    TrialFunction,
    )

#Local libraries and methods
from local_library1 import (
    class1,
    method1, 
    method2,
    )
import local_library2 as lib2
```
## Classes
The code is to be written object oriented. The names of private attributes and functions are preceeded by a double underscore '\_\_', the ones of protected attributes and function by a single underscore '\_'. By default, all attributes of a class are private. If public access needs to be given to a private variable, this is handled via the property decorator:

**CORRECT:**
```python
    #In the class
	@property
	def privateVariable(self):
		return self.__privateVariable
		
    #Access to the attribute from a class-object
    theObjectsPrivateVariable = ObjectOfClass.privateVariable
```
In the body of the class, first private, then protected and then public functions are defined. Within these groups, the functions are ordered alphabetically.

## Documentation \& Comments
### Comments
Comments are meant to clarify the code, where a future reader may have trouble understanding it. Nevertheless, it is important to not clutter the code with comments, which are redundant to the code, or confusing. Please follow the following rules when commenting the code (based on [Ellen Spertus](https://stackoverflow.blog/2021/12/23/best-practices-for-writing-code-comments/)):

1. Comments should not duplicate the code
2. Avoid unnecessary comments by writing clean code using descriptive variable and function names. It is the worst to write bad code without comments, it is better (still bad) to write bad code with comments, and it is best to write clear code, which makes comments unnecessary.
3. Comment code which might be misleadingly considered unneeded or redundant by a reader.
4. Add comments, where the code should be adapted. Start these comments with 'TODO'.

**WRONG:**
```python
	# Add one to i
    i += 1;         
```

### Docstrings
Docstrings are to be added for classes and functions according to the [numpy docstring style](https://numpydoc.readthedocs.io/en/latest/format.html).
Examples can be found [here](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).



