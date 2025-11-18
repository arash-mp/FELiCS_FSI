# How to import and export a h5 file

Before you start, you should install all necessary FELiCS packages, and install FELiCS as a package itself. The [installation guide](https://felics2-0-laboratory-for-flow-instabilities-and--112f91add91c56.gitlab-pages.tu-berlin.de/installation_guide.html) explains those steps in detail.

After installing FELiCS you can start writing your first FELiCS script. 

Here, we have used gmsh to create a square mesh, which we exported in ASCII format. If you need a more detailed explanation you can have a look [here](https://felics2-0-laboratory-for-flow-instabilities-and--112f91add91c56.gitlab-pages.tu-berlin.de/Tutorials/cylinder_wake.html).

**Step 1: Create a FELiCSMesh**

The FELiCSMesh is a wrapper for a dolfinx mesh, that contains additional attributes and methods, e.g. the coordinate system of your case. 


```python
from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh

meshFileName = "./square_mesh.msh"
coordinateSystemName = "Cartesian"

# Create a FELiCS mesh from the gmsh file
mesh = FELiCSMesh(coordinateSystemName=coordinateSystemName, meshFileName=meshFileName)
```

    [38;20mInfo     | FELiCSMesh.py          | __init__                   (line 78  ) : Opening mesh file: ./square_mesh.msh[0m
    [38;20mInfo     | FELiCSMesh.py          | __init__                   (line 86  ) : Mesh contains 513 nodes and 1024 elements[0m


**Step 2: Create a scalar function space**

Next we create a dolfinx function space from our FELiCSMesh object with the polynomial order and the dimension of the space. In FELiCS, only continuous Lagrange elements are used. Because our space should be scalar, we choose dim=1.


```python
from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace

space = getFELiCSSpace(mesh, order=2, dim=1)
```

**Step 3: Create fields and import from H5 files**

Now things are getting interesting. The Field object can store fields on our function space, as the name suggests. It contains a dolfinx function, as well as a FELiCS tensor object (which we use to make our expressions coordinatesystem independent). We can import data to our Field from a h5 file, but also export our created data to h5. Fields can be naturally added to eachother using the standard python operators. Make sure, that the fields are defined on the same space, otherwise the sum cannot be computed.


```python
from FELiCS.Fields.Field import Field

Field1 = Field(space, mesh=mesh)
Field1.importH5("function_values_quadratic")
Field2 = Field(space, mesh=mesh)
Field2.importH5("function_values_sine")

Field3 = Field1 + Field2

Field3.exportH5("field_sum")

```

**Visualize the results**

The resulting "xdmf" file can be load e.g. into Paraview. Our field now looks like this: (TODO: get plot from Field)
