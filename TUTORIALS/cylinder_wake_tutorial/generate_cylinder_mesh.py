#!/usr/bin/env python3
"""
Mesh Generator for Flow Around Sphere Wake analysis
Author      : [Simon Demange]
Date        : [15.05.2025]
Modified    : [Anant Talasikar, 21.05.2025]
Description:

This script generates a 2D mesh for the linear analysis of the flow around a sphere 
(represented as a circle in 2D), with refined mesh regions around and behind the sphere 
to capture wake effects. The mesh is designed for use with FELiCS.

The geometry consists of:
1. A large rectangular domain with specific dimensions
2. A circle representing a sphere cross-section
3. Two refinement regions with progressively finer mesh sizing:
    - A larger wake region behind the sphere
    - A finer region immediately surrounding and behind the sphere

GMSH methods used:
- gmsh.initialize(): Initializes the gmsh API
- gmsh.model.geo.addPoint(x, y, z, meshSize, tag):                              Adds a point at coordinates (x,y,z) with specified mesh size
- gmsh.model.geo.addLine(startPointTag, endPointTag, tag):                      Creates a line between two points
- gmsh.model.geo.addCircleArc(startPointTag, centerPointTag, endPointTag, tag): Creates a circular arc
- gmsh.model.geo.synchronize():                                                 Synchronizes the CAD model with the gmsh model
- gmsh.model.addPhysicalGroup(dim, tags, tag, name):                            Creates named physical groups for boundary conditions
- gmsh.model.geo.addCurveLoop(curveTags, tag):                                  Creates a closed loop from curves
- gmsh.model.geo.addPlaneSurface(wireTags, tag):                                Creates a surface from closed-loop boundaries
- gmsh.model.mesh.generate(dim):                                                Generates a mesh with specified dimension
- gmsh.option.setNumber(name, value):                                           Sets gmsh options
- gmsh.write(filename):                                                         Writes the mesh to a file
- gmsh.finalize():                                                              Finalizes and cleans up the gmsh API

Output: sphere_wake.msh - A 2D mesh file for FELiCS, in msh VERSION 2.0 ASCII format
"""
class create_custom_mesh:
    
    def __init__(self, refmiX, refmaX, refmaY, refmi2X, refma2X, refma2Y, grf, rfc, rff, rffest):
        self.RefineMinX                 = refmiX
        self.RefineMaxX                 = refmaX
        self.RefineMaxY                 = refmaY
        self.Refine2MinX                = refmi2X
        self.Refine2MaxX                = refma2X
        self.Refine2MaxY                = refma2Y 

        # Mesh sizes
        self.globalRefinementFactor     = grf
        self.RefinementFactorCorarse    = rfc
        self.RefinementFactorFine       = rff
        self.RefinementFactorFinest     = rffest
        # dxCoarse                    = RefinementFactorCorarse*globalRefinementFactor
        # dxFine                      = RefinementFactorFine*globalRefinementFactor
        # dxFinest                    = RefinementFactorFinest*globalRefinementFactor



    def gen_cylinder_mesh(MeshDef):
        # Package import
        import gmsh
        
        # General geometry parameters
        MinX        = -100
        MaxX        = 200
        MinY        = 0
        MaxY        = 25
        R           = 0.5
        dxCoarse    = MeshDef.RefinementFactorCorarse* MeshDef.globalRefinementFactor
        dxFine      = MeshDef.RefinementFactorFine*MeshDef.globalRefinementFactor
        dxFinest    = MeshDef.RefinementFactorFinest*MeshDef.globalRefinementFactor 


        # Add the points to model
        # Format is            (x,                      y,                      z,  mesh size,  index)
        gmsh.initialize()
        gmsh.model.geo.addPoint(MinX,                   MinY,                   0,  dxCoarse,    1)
        gmsh.model.geo.addPoint(MeshDef.RefineMinX,     MinY,                   0,  dxCoarse,    2)
        gmsh.model.geo.addPoint(MeshDef.RefineMinX,     MeshDef.RefineMaxY,     0,  dxCoarse,    3)
        gmsh.model.geo.addPoint(MeshDef.RefineMaxX,     MeshDef.RefineMaxY,     0,  dxCoarse,    4)
        gmsh.model.geo.addPoint(MeshDef.RefineMaxX,     MinY,                   0,  dxCoarse,    5)
        gmsh.model.geo.addPoint(MaxX,                   MinY,                   0,  dxCoarse,    6)
        gmsh.model.geo.addPoint(MaxX,                   MaxY,                   0,  dxCoarse,    7)
        gmsh.model.geo.addPoint(MinX,                   MaxY,                   0,  dxCoarse,    8)

        gmsh.model.geo.addPoint(-R,                     MinY,                   0,  dxFinest,    9)
        gmsh.model.geo.addPoint(0,                      R,                      0,  dxFinest,    10)
        gmsh.model.geo.addPoint(R,                      MinY,                   0,  dxFinest,    11)
        gmsh.model.geo.addPoint(0,                      0,                      0,  dxFinest,    12)

        gmsh.model.geo.addPoint(MeshDef.Refine2MinX,    MinY,                   0,  dxFine,      13)
        gmsh.model.geo.addPoint(MeshDef.Refine2MinX,    MeshDef.Refine2MaxY,    0,  dxFine,      14)
        gmsh.model.geo.addPoint(MeshDef.Refine2MaxX,    MeshDef.Refine2MaxY,    0,  dxFine,      15)
        gmsh.model.geo.addPoint(MeshDef.Refine2MaxX,    MinY,                   0,  dxFine,      16)

        # Define the lines
        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.addLine(2, 3, 2)
        gmsh.model.geo.addLine(3, 4, 3)
        gmsh.model.geo.addLine(4, 5, 4)
        gmsh.model.geo.addLine(5, 6, 5)
        gmsh.model.geo.addLine(6, 7, 6)
        gmsh.model.geo.addLine(7, 8, 7)
        gmsh.model.geo.addLine(8, 1, 8)
        gmsh.model.geo.addLine(13, 9, 9)

        # Define the are-circles for the cylinder
        gmsh.model.geo.addCircleArc(9, 12, 10, 10)
        gmsh.model.geo.addCircleArc(10, 12, 11, 11)

        # More lines 
        gmsh.model.geo.addLine(11, 16, 12)
        gmsh.model.geo.addLine(16, 15, 13)
        gmsh.model.geo.addLine(15, 14, 14)
        gmsh.model.geo.addLine(14, 13, 15)
        gmsh.model.geo.addLine(2, 13, 16)
        gmsh.model.geo.addLine(16, 5, 17)

        # Define the physical groups that will be used for 
        # boundary conditions in FELiCS
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(1, [8],                 1, "Inlet")
        gmsh.model.addPhysicalGroup(1, [1,16,9,12,17,5],    2, "Symmetry")
        gmsh.model.addPhysicalGroup(1, [6],                 3, "Outlet")
        gmsh.model.addPhysicalGroup(1, [7],                 4, "Top")
        gmsh.model.addPhysicalGroup(1, [10,11],             5, "Wall")

        # Add the closed contours and surface
        gmsh.model.geo.addCurveLoop([1,2,3,4,5,6,7,8], 100)
        gmsh.model.geo.addPlaneSurface([100], 101)
        gmsh.model.geo.addCurveLoop([9,10,11,12,13,14,15], 200)
        gmsh.model.geo.addPlaneSurface([200], 201)
        gmsh.model.geo.addCurveLoop([16,-15,-14,-13,17,-4,-3,-2], 300)
        gmsh.model.geo.addPlaneSurface([300], 301)

        # Merge all surfaces into one
        gmsh.model.addPhysicalGroup(2, [101,201,301], 6, "all") # domain

        # Finalize geometry and mesh
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)

        # Define type of export file
        # gmsh.model.setColor(2, 1, 0, 0, 0)
        # gmsh.option.setNumber("Mesh.ElementEdges", 1)  # Make sure edges are shown
        # gmsh.option.setNumber("Mesh.ElementEdgesColor", 1)  # Set edge color to black (index 1)
        # gmsh.option.setNumber("Mesh.ColorCarousel", 1)  # Set edge color to black (index 1)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write("cylinder_wake.msh")

        gmsh.finalize()

