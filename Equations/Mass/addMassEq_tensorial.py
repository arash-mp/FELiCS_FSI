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
from functions import printWarning, printError

def addMassEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    mass conservation equation, in tensorial framework.
    
    NOTE:   The "FLAG_TENS" is used to activate the tensor
            framework in specific terms of the eq. It is possible
            to combine tensorial terms with index-notation ones
    '''

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')        

    # ------------------------ Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED FOR CARTESIAN COORDINATES, should be moved to param
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        rhou_fluc_tens = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], 0.0)), coord_sys)
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # HARDCODED FOR CYL COORD IN FELiCS, DEFINED AS (Z, R, PHI)!
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        rhou_fluc_tens = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], fluc.rhou[2])), coord_sys)
        
    else:
        printWarning('Coord. syst not yet implemented in tensor framework.')

    rho_fluc_tens = Tensor(fluc.rho, coord_sys)
    x_tens = Tensor(X, coord_sys)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)

    # ------------------------ Time derivative term used
    # Only if density fluctuations are considered
    FLAG_TENS = True
    
    if 'rho' in param.Case.getTransportedQuantityList():
        if FLAG_TENS:
            printWarning('Density fluctuation term not yet validated in tensor framework.')
            self.B_vf.add( (rho_fluc_tens*iConj(x_tens)).ufl_tens*coord_sys.J_hat*dx)
            
        else:
            # -- > Previous implementation
            self.B_vf.add( fluc.rho*self.R*conj(X)*dx)


    # ------------------------ Advection terms
    # The advection term is integrated by parts
    FLAG_TENS = True
    
    if FLAG_TENS:
        # -- > Tensor implementation
        # Volume term from IbP
        self.A_vf.add( (1j*iDot(iGrad(iConj(x_tens)), rhou_fluc_tens)).ufl_tens*coord_sys.J_hat*dx)
        
        # Boundary term from IbP
        self.A_vf.add( (-1j*iDot(nbc_tens, rhou_fluc_tens*iConj(x_tens)) ).ufl_tens*coord_sys.J_hat*self.all_ds)
        
        # Extra term for m > 0
        if (param.Case.CoordinateSystem in ['Cylindrical']) and (not param.Case.m == 0):
            printWarning('m > 0 term not yet implemented in tensor framework. Results are wrong!')
            # We can't just add the terms for m > 1 in the tensor framework 
            self.A_vf.add(conj(X)*mean.rho*self.m*fluc.u[2]*dx)
            self.A_vf.add(conj(X)*fluc.rho*self.m*mean.u[2]*dx)
        
    else:
        # -- > Previous implementation
        I = Identity(fluc.u.geometric_dimension())
        
        # Volume term from IbP
        self.A_vf.add(1j * conj(X).dx(i)*self.R*fluc.rhou[j]*I[i,j]*dx)
        
        # Boundary term from IbP
        self.A_vf.add(1j * -self.R*fluc.rhou[j]*self.n_BC[i]*conj(X)*I[i,j]*self.all_ds)
        
        # Extra term for m > 0 
        if (param.Case.CoordinateSystem in ['Cylindrical']) and (not param.Case.m == 0):
            self.A_vf.add(conj(X)*mean.rho*self.m*fluc.u[2]*dx)
            self.A_vf.add(conj(X)*fluc.rho*self.m*mean.u[2]*dx)

    # ------------------------ BC term for Input/Output analysis
    FLAG_TENS = True
    
    if param.Case.AnalysisMode in ['Input-Output']:
        if FLAG_TENS:
            printWarning('Input-Output analysis not implemented in tensor framework. Currently relies on index notation.')

        # -- > Previous implementation
        # Iterate through all boundaries, at which forcing is applied
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            # First subtract the part added in a few lines above...
            self.A_vf.add(1j * self.R*inner(fluc.u,self.n_BC)*conj(X)*self.ds(boundary_index))
            # Then add the forcing at the inlet...
            self.A_vf.add(-self.R*inner(mean.u_forcing_r,self.n_BC)\
                        *conj(X)*self.ds(boundary_index))
            self.A_vf.add(1j * -self.R*inner(mean.u_forcing_i,self.n_BC)\
                        *conj(X)*self.ds(boundary_index))
