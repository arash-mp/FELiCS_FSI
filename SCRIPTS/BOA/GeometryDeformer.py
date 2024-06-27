from abc import ABC, abstractmethod
import dolfinx
import numpy as np
import copy
from FELiCS.Misc.functions import printWarning


class GeometryDeformer(ABC):

    def __init__(self, mesh, facets):
        self.mesh       = mesh
        self.facets     = facets
        self.isDeformed = False

        self.x         = mesh.geometry.x
        self.x_dolfinx = mesh.dolfinxMesh.geometry.x
        self.x_cpp     = mesh._cpp_object.geometry.x

        bc_facet_dofs  = facets.indices[facets.values==1001]
        self.bc_dofs   = dolfinx.mesh.compute_incident_entities(self.mesh.topology,bc_facet_dofs, self.mesh.topology.dim-1, 0)

        
        self.x_save   = copy.deepcopy(self.x[self.bc_dofs])
        #self.x_save = [None]*len(self.x[self.bc_dofs])
        #i=0
        #for dof in self.bc_dofs:
        #    self.x_save[i]    = np.empty(3)
        #    self.x_save[i]    = copy.deepcopy(self.x[dof])
        #    i+=1

    @abstractmethod
    def getNumberOfParameters(self):
        pass

    @abstractmethod
    def deformMesh(self, parameters):
        pass

    def restoreMesh(self):
        i=0
        for dof in self.bc_dofs:
            self.x[dof]         = self.x_save[i]
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
        from   geomdl import utilities

        if self.isDeformed:
            #TODO: printWarning
            printWarning('CAUTION: deforming an already deformed mesh...')

        curve = BSpline.Curve()
        
        # Set degree
        curve.degree = 3
        
        # Set control points
        curve.ctrlpts = [[-0.5, 0, 0], [-0.25, parameters[0], 0], [0.25, parameters[1], 0], [0.5, 0, 0]]
        
        # Set knot vector
        curve.knotvector = [0, 0, 0, 0, 1, 1, 1, 1]
        
        curve.delta = 1.e-4
        curve_points = curve.evalpts
        N = len(curve_points)
        x_val = np.empty(N)
        y_val = np.empty(N)
        for i in range(N):
            x_val[i] = curve_points[i][0]
            y_val[i] = curve_points[i][1]
        
        spline = interpolate.CubicSpline(x_val, y_val)
        
        for dof in self.bc_dofs:
            self.x[dof][1]         = spline(self.x[dof][0])
            self.x_dolfinx[dof][1] = spline(self.x[dof][0])
            self.x_cpp[dof][1]     = spline(self.x[dof][0])

        self.isDeformed=True


class ProfileBSpline8Pts(GeometryDeformer):

    def __init__(self,mesh,facets):
        super().__init__(mesh, facets)
        self.NumberOfParams = 8

    def getNumberOfParameters(self):
        return self.NumberOfParams

    def deformMesh(self, parameters):
        import scipy.interpolate as interpolate
        import geomdl as g
        from   geomdl import BSpline
        from   geomdl import utilities

        if self.isDeformed:
            #TODO: printWarning
            printWarning('CAUTION: deforming an already deformed mesh...')

        curve = BSpline.Curve()
        
        # Set degree
        curve.degree = 3
        
        # Set control points
#        curve.ctrlpts = [[0.10000000000000000555, -0.00499999999999998536, 0],\
#                         [0.09777777777777778290, -0.00242300886715957411, 0],\
#                         [0.09555555555555556024, -0.00042898338465717500, 0],\
#                         [0.09333333333333333759,  0.00114762923408253878, 0],\
#                         [0.09111111111111111494,  0.00239516041194041962, 0],\
#                         [0.08888888888888890616,  0.00336637871645898909, 0],\
#                         [0.08666666666666668351,  0.00409472049133493068, 0],\
#                         [0.08444444444444446085,  0.00460176647276058781, 0],\
#                         [0.08222222222222223820,  0.00490103870112774259, 0],\
#                         [0.08000000000000000167,  0.00500000000000000097, 0]]
#                

        ## only y: (NOTE: because for the CubicSpline x has to be increasing, I changed the order of the points, which may look weird. The parameters do have the same order as the points in the geo file.)
        curve.ctrlpts = [[0.08000000000000000167,  0.00500000000000000097, 0],\
                         [0.08222222222222223820, parameters[7], 0],\
                         [0.08444444444444446085, parameters[6], 0],\
                         [0.08666666666666668351, parameters[5], 0],\
                         [0.08888888888888890616, parameters[4], 0],\
                         [0.09111111111111111494, parameters[3], 0],\
                         [0.09333333333333333759, parameters[2], 0],\
                         [0.09555555555555556024, parameters[1], 0],\
                         [0.09777777777777778290, parameters[0], 0],\
                         [0.10000000000000000555, -0.00499999999999998536, 0]]

        # Set knot vector
        curve.knotvector = utilities.generate_knot_vector(curve.degree, curve.ctrlpts_size) 

        curve.delta = 1.e-4
        curve_points = curve.evalpts
        N = len(curve_points)
        x_val = np.empty(N)
        y_val = np.empty(N)
        for i in range(N):
            x_val[i] = curve_points[i][0]
            y_val[i] = curve_points[i][1]
        
        spline = interpolate.CubicSpline(x_val, y_val)
        
        for dof in self.bc_dofs:
            self.x[dof][1]         = spline(self.x[dof][0])
            self.x_dolfinx[dof][1] = spline(self.x[dof][0])
            self.x_cpp[dof][1]     = spline(self.x[dof][0])

        self.isDeformed=True




