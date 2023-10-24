from ufl import (
    inner,
    dx,
    Identity,
    i,
    j,
    conj,
)
from tensor_utils import (
    Tensor,
    CoordinateSystem,
    as_vector,
    iGrad,
    iDot,
    iConj
)
from functions import printWarning

def addMassEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    mass conservation equation, in tensorial framework.
    '''

    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printWarning('Discontinuous Galerkin not implemented in tensorial framework.')

    if param.Case.CoordinateSystem =='Cylindrical':
        printWarning('Cyl. coord. not yet tested in tensor framework.')

    if not param.Case.m == 0:
        printWarning('m > 0 not yet implemented in tensor framework. \
                     For cyl. coord. it relies on previous ufl implementation.')

    if param.Case.AnalysisMode in ['Input-Output']:
        printWarning('Input-Outpout terms not yet implemented in tensor framework. \
                     It relies on previous ufl implementation.')

    if 'rho' in param.Case.getTransportedQuantityList():
        printWarning('Density fluctuations not yet tested in tensor framework.')

    # ------------------------ Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED FOR CARTESIAN COORDINATES, should be moved to param
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        rho_fluc_tens = Tensor(fluc.rho, coord_sys)
        rhou_fluc_tens = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], 0.0)), coord_sys)
        x_tens = Tensor(X, coord_sys)
        nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # HARDCODED FOR CYL COORD IN FELiCS, DEFINED AS (Z, R, PHI)!
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        rho_fluc_tens = Tensor(fluc.rho, coord_sys)
        rhou_fluc_tens = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], fluc.rhou[2])), coord_sys) # We should consider the three velocity components but this seem to add some error!
        x_tens = Tensor(X, coord_sys)
        nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)
        
    else:
        printWarning('Coord. syst not yet implemented in tensor framework.')


    # ------------------------ Time derivative term used
    # (if density fluctuations are considered)
    # THIS TERM WAS NOT TESTED YET!
    if 'rho' in param.Case.getTransportedQuantityList():
        # -- > Previous implementation
        #self.B_vf.add( fluc.rho*self.R*conj(X)*dx)
        # -- > Tensor implementation (not working for cyl coord)
        self.B_vf.add( (rho_fluc_tens*iConj(x_tens)).ufl_tens*coord_sys.J_hat*dx)


    # ------------------------ Advection terms
    # -- > Previous implementation
    I = Identity( fluc.u.geometric_dimension() )
    #self.A_vf.add(1j * conj(X).dx(i)*coord_sys.J_hat*fluc.rhou[j]*I[i,j]*dx)
    # -- > Tensor implementation (not working for cyl coord)
    self.A_vf.add( (1j*iDot(iGrad(iConj(x_tens)), rhou_fluc_tens)).ufl_tens*coord_sys.J_hat*dx)

    # An integration by part is used to make BC term appear?
    # -- > Previous implementation
    #self.A_vf.add(1j * -self.R*fluc.rhou[j]*self.n_BC[i]*conj(X)*I[i,j]*self.all_ds)
    # -- > Tensor implementation (not working for cyl coord)
    self.A_vf.add( (-1j*iDot(nbc_tens, rhou_fluc_tens*iConj(x_tens)) ).ufl_tens*coord_sys.J_hat*self.all_ds)

    # m > 0 NOT IMPLEMENTED IN TENSOR FRAMEWORK YET
    if (param.Case.CoordinateSystem in ['Cylindrical']) and (not param.Case.m == 0):
        self.A_vf.add(conj(X)*mean.rho*self.m*fluc.u[2]*dx)
        self.A_vf.add(conj(X)*fluc.rho*self.m*mean.u[2]*dx)

    # In case of Input-Output analysis, we must adapt the boundary terms...
    # NOT YET TESTED FOR TENSOR FRAMEWORK
    if param.Case.AnalysisMode in ['Input-Output']:
        # Iterate through all boundaries, at which forcing is applied
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            # First subtract the part added in a few lines above...
            self.A_vf.add(1j * self.R*inner(fluc.u,self.n_BC)*conj(X)*self.ds(boundary_index))
            # Then add the forcing at the inlet...
            self.A_vf.add(-self.R*inner(mean.u_forcing_r,self.n_BC)\
                          *conj(X)*self.ds(boundary_index))
            self.A_vf.add(1j * -self.R*inner(mean.u_forcing_i,self.n_BC)\
                          *conj(X)*self.ds(boundary_index))
