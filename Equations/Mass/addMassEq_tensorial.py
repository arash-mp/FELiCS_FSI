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
        
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(x)
    '''

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')        


    # ------------------------ Define the tensorial operators
    # Coordinate system
    coord = self.coord_sys
    
    # Vector quantities
    if fluc.rhou.ufl_shape[0] == 2:
        rhou_f = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], 0.0)), coord)
        if param.Case.AnalysisMode in ['Input-Output']:
            forcing_comp = mean.u_forcing_r + 1j*mean.u_forcing_i
            forcing_tens = Tensor(as_vector((forcing_comp[0], forcing_comp[1], 0.0)), coord)
            
    elif fluc.rhou.ufl_shape[0] == 3:
        rhou_f = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], fluc.rhou[2])), coord)
        if param.Case.AnalysisMode in ['Input-Output']:
            forcing_comp = mean.u_forcing_r - 1j*mean.u_forcing_i
            forcing_tens = Tensor(as_vector((forcing_comp[0], forcing_comp[1], forcing_comp[2])), coord)
            
    else:
        printError('--> Mass eq: rhou has neither 2 or 3 dimensions: not implemented.')
    
    # Scalar quantities    
    rho_f = Tensor(fluc.rho, coord)
    x_tens = Tensor(X, coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)


    # ------------------------ Time derivative term used
    # Only if density fluctuations are considered
    FLAG_TENS = True
    if 'rho' in param.Case.getTransportedQuantityList():
        if FLAG_TENS:
            self.B_vf.add( (rho_f*iConj(x_tens)).ufl_tens*coord.J_hat*dx)
            
        else:
            # -- > Previous implementation
            self.B_vf.add( fluc.rho*self.R*conj(X)*dx)


    # ------------------------ Advection terms
    # The advection term is integrated by parts
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation
        # Volume term from IbP
        self.A_vf.add(( 1j*iDot(iGrad(iConj(x_tens),-self.m),rhou_f) ).ufl_tens*coord.J_hat*dx)
        # Boundary term from IbP
        self.A_vf.add(( -1j*iDot(nbc_tens,rhou_f*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
        
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
            # Iterate through all boundaries, at which forcing is applied
            for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
                # First subtract the boundary term from advection
                # Correction wrt to index notation: We need to remove the rho*u term, not just the u! 
                self.A_vf.add(( 1j*iDot(nbc_tens,rhou_f*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(boundary_index))
                # Then add the forcing at the boundary
                self.A_vf.add(( -1*iDot(nbc_tens,forcing_tens)*iConj(x_tens) ).ufl_tens*coord.J_hat*self.ds(boundary_index))
                
        else:
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
