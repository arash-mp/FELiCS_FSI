```py
import gmsh
```
After importing the python-gmsh package the following methods become accessible:
- ```gmsh.initialize()```: Initializes the gmsh API
- ```gmsh.model.geo.addPoint(x, y, z, meshSize, tag)```:                              Adds a point at coordinates (x,y,z) with specified mesh size
- ```gmsh.model.geo.addLine(startPointTag, endPointTag, tag)```:                      Creates a line between two points
- ```gmsh.model.geo.addCircleArc(startPointTag, centerPointTag, endPointTag, tag)```: Creates a circular arc
- ```gmsh.model.geo.synchronize()```:                                                 Synchronizes the CAD model with the gmsh model
- ```gmsh.model.addPhysicalGroup(dim, tags, tag, name)```:                            Creates named physical groups for boundary conditions
- ```gmsh.model.geo.addCurveLoop(curveTags, tag)```:                                  Creates a closed loop from curves
- ```gmsh.model.geo.addPlaneSurface(wireTags, tag)```:                                Creates a surface from closed-loop boundaries
- ```gmsh.model.mesh.generate(dim)```:                                                Generates a mesh with specified dimension
- ```gmsh.option.setNumber(name, value)```:                                           Sets gmsh options
- ```gmsh.write(filename)```:                                                         Writes the mesh to a file
- ```gmsh.finalize()```:                                                              Finalizes and cleans up the gmsh API

Output: sphere_wake.msh - A 2D mesh file for FELiCS, in msh VERSION 2.0 ASCII format
