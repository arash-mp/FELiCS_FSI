# Calculate the L2-norm of a field

**Step 1: Create FELiCSMesh, FELiCSSpace and Field**



```python
from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh

meshFileName = "./square_mesh.msh"
coordinateSystemName = "Cartesian"

# Create a FELiCS mesh from the gmsh file
mesh = FELiCSMesh(coordinateSystemName=coordinateSystemName, meshFileName=meshFileName)

from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace

scalarSpace = getFELiCSSpace(mesh, order=2, dim=1)

from FELiCS.Fields.Field import Field

phi = Field(scalarSpace, mesh=mesh)
phi.importH5("function_values_sine")

```

**Step 2: Calculate the L2-norm of a field**

Now we will calculate an integral quantity of our field. We will calculate the L2-norm of our field over the whole domain.

The L2-norm is defined as:
$$||\phi||_{L^2} = \sqrt{\int_\Omega |\phi|^2 \, d\Omega}$$

For this, the `Field` class has the method `calculateL2Norm()` which computes this integral over the entire mesh domain.



```python
u_norm = phi.calculateL2Norm()

print(f"L2-norm of the field: {u_norm}")
```
