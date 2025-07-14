# How to: First steps in FELiCS 

[Link to installation guide]
How to import FELiCS...


```python
The <a href="https://git.tu-berlin.de/laboratory-for-flow-instabilities-and-dynamics/felics2.0/-/blob/main/DOCUMENTATION/installation_guide.md"> installation guide </a> explains how to create the conda environment, download FELiCS and install FELiCS as a package.

After installing FELiCS you can start writing your first FELiCS script.
```

We use gmsh to create a mesh. [TODO: link to tutorial]
First we create a FELiCSMesh object from a gmsh mesh file. The FELiCSMesh is a wrapper for a dolfinx mesh, that contains additional attributes and methods, e.g. the coordinate system of your case.



```python
from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh

meshFileName = "./square_mesh.msh"
coordinateSystem = "Cartesian"

# Create a FELiCS mesh from the gmsh file
mesh = FELiCSMesh(coordinateSystem=coordinateSystem, meshFileName=meshFileName)
```

Next we create a FEM function space from our FELiCSMesh object.


```python
from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace

space = getFELiCSSpace(mesh, order=2, dim=1)


```

Now things are getting interesting. The Field object can store fields on our function space, as the name suggests. We can import data to our Field from a h5 file, but also export our created data to h5. Fields can be naturally added to eachother using the standard python operators. Make sure, that the fields are defined on the same space, otherwise the fields cannot be added together.


```python
from FELiCS.Fields.Field import Field

Field1 = Field(space, mesh=mesh)
Field1.importH5("function_values_quadratic")
Field2 = Field(space, mesh=mesh)
Field2.importH5("function_values_sine")

Field3 = Field1 + Field2

Field3.exportH5("field_sum")

```
