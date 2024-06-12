from abc import ABC, abstractmethod
import dolfinx
import numpy as np

class GeometryDeformer(ABC):

    def __init__(self, mesh, facets):
        self.mesh       = mesh
        self.facets     = facets
        self.isDeformed = False

        self.x = mesh.geometry.x
        bc_facet_dofs = facets.indices[facets.values==1001]
        self.bc_dofs  = dolfinx.mesh.compute_incident_entities(self.mesh.topology,bc_facet_dofs, self.mesh.topology.dim-1, 0)
        
        self.x_save = [None]*len(self.x[self.bc_dofs])
        i=0
        for dof in self.bc_dofs:
            self.x_save[i] = np.empty(3)
            self.x_save[i][:] = self.x[dof][:]
            i+=1

    @abstractmethod
    def getNumberOfParameters(self):
        pass

    @abstractmethod
    def deformMesh(self, parameters):
        pass

    def restoreMesh(self):
        i=0
        for dof in self.bc_dofs:
            self.x[dof][:] = self.x_save[i][:]
            i+=1
        self.isDeformed = False




class CylinderBSpline2Pts(GeometryDeformer):

    def __init__(self,mesh,facets):
        super().__init__(mesh, facets)
        self.NumberOfParams = 2

    def getNumberOfParameters(self):
        return self.NumberOfParams

    def deformMesh(self, parameters):
        import scipy.interpolate as interpolate
        import geomdl as g
        from   geomdl import BSpline

        curve = BSpline.Curve()
        
        # Set degree
        curve.degree = 3
        
        # Set control points
        curve.ctrlpts = [[-0.5, 0, 0], [-0.25, parameters[0], 0], [0.25, parameters[1], 0], [0.5, 0, 0]]
        
        # Set knot vector
        curve.knotvector = [0, 0, 0, 0, 1, 1, 1, 1]
        
        curve.delta = 1.e-3
        curve_points = curve.evalpts
        N = len(curve_points)
        x_val = np.empty(N)
        y_val = np.empty(N)
        for i in range(N):
            x_val[i] = curve_points[i][0]
            y_val[i] = curve_points[i][1]
        
        spline = interpolate.CubicSpline(x_val, y_val)
        
        for dof in self.bc_dofs:
            self.x[dof][1] = spline(self.x[dof][0])




