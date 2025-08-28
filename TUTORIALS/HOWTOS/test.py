from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh

meshFileName = "./square_mesh.msh"

coordinateSystemName = "Cartesian"

# Create a FELiCS mesh from the gmsh file
mesh = FELiCSMesh(coordinateSystemName=coordinateSystemName, meshFileName=meshFileName)

from FELiCS.SpaceDisc.FEMSpaces import getFELiCSSpace

scalarSpace = getFELiCSSpace(mesh, order=2, dim=1)
vectorSpace = getFELiCSSpace(mesh, order=2, dim=2)

from FELiCS.Fields.Field import Field

phi = Field(vectorSpace, mesh=mesh)
phi.importH5("function_values_2d")
vorticitySol = Field(scalarSpace, mesh=mesh)
vorticitySol.importH5("vorticity_solution")


import numpy as np

vorticity = phi.getVorticityField()

# print(cg[0].getListOfSingleFields())
# print(cg[0].getListOfSingleFields()[1].getCoefficientArray())
# print(cg[1].getListOfSingleFields()[0].getCoefficientArray())


print(vorticity.getCoefficientArray())
print(vorticitySol.getCoefficientArray())
print(np.linalg.norm(vorticity.getCoefficientArray() - vorticitySol.getCoefficientArray()))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.tri as tri

# Calculate the magnitude of the vector field
function_values_x = phi.getListOfSingleFields()[0].getCoefficientArray()
function_values_y = phi.getListOfSingleFields()[1].getCoefficientArray()
magnitude = np.sqrt(function_values_x**2 + function_values_y**2)

# Extract x and y coordinates from dof_coordinates
dof_coordinates = vectorSpace.tabulate_dof_coordinates()
x_coords = dof_coordinates[:, 0]
y_coords = dof_coordinates[:, 1]

vorticity_values = vorticity.getCoefficientArray()
# Create triangulation for irregular grid
triang = Triangulation(x_coords, y_coords)

# Create the contour plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Filled contour plot
contour_filled = ax1.tricontourf(triang, magnitude, levels=20, cmap='viridis')
ax1.set_title('Magnitude of Vector Field (Filled Contours)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_aspect('equal')
plt.colorbar(contour_filled, ax=ax1, label='Magnitude')

print(magnitude)
print(vorticity_values)
print(np.max(vorticity_values))
print(np.min(vorticity_values))

# Plot 1: Filled contour plot
vorticity_filled = ax2.tricontourf(triang, vorticity_values, levels=20, cmap='viridis')
ax2.set_title('Vorticity')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_aspect('equal')
plt.colorbar(vorticity_filled, ax=ax2, label='Vorticity')

plt.tight_layout()
plt.show()