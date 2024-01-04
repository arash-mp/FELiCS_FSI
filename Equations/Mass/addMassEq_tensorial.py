from ufl import dx
from tensorUtils import (
    Tensor,
    as_vector,
    iGrad,
    iDot,
    iConj
)
from functions import printDebug, printError

def addMassEq(self,fluc,X_in,mean,param):
    '''
    This function builds the weak form of the linearized
    mass conservation equation, in tensorial framework.
        
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
    
    ## Vector quantities
    #if fluc.rhou.ufl_shape[0] == 2:
    #    printDebug(param.debug,'--> Mass eq: fluctuations and mean flow 2D')
    #    #fluc_rhou = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], 0.0)), coord)
    #    fluc_rhou = Tensor(fluc.rhou, coord)
    #    #input(len(fluc.p.ufl_shape))
    #    #input(len(fluc.rhou.ufl_shape))
    #    if param.Case.AnalysisMode in ['Input-Output']:
    #        forcing_comp = mean.u_forcing_r + 1j*mean.u_forcing_i
    #        forcing_tens = Tensor(as_vector((forcing_comp[0], forcing_comp[1], 0.0)), coord)
    #        
    #elif fluc.rhou.ufl_shape[0] == 3:
    #    printDebug(param.debug,'--> Mass eq: fluctuations and mean flow 3D')
    #    fluc_rhou = Tensor(as_vector((fluc.rhou[0], fluc.rhou[1], fluc.rhou[2])), coord)
    #    if param.Case.AnalysisMode in ['Input-Output']:
    #        forcing_comp = mean.u_forcing_r + 1j*mean.u_forcing_i
    #        forcing_tens = Tensor(as_vector((forcing_comp[0], forcing_comp[1], forcing_comp[2])), coord)
    #        
    #else:
    #    printError('--> Mass eq: rhou has neither 2 or 3 dimensions: not implemented.')
    fluc_rhou = fluc.rhou 
    # Scalar quantities    
    fluc_rho = fluc.rho
    #X = Tensor(X_in, coord,containsTestFunction = True)
    X = X_in
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)


    # ------------------------ Time derivative term used
    # Only if density fluctuations are considered
    if 'rho' in param.Case.getTransportedQuantityList():
        self.B_vf.add( (fluc_rho*iConj(X)).ufl_tens*coord.J_hat*dx)


    # ------------------------ Advection terms
    # The advection term is integrated by parts
    # Volume term from IbP
    self.A_vf.add(( 1j*iDot(iGrad(iConj(X), self.m),fluc_rhou) ).ufl_tens*coord.J_hat*dx)
    # Boundary term from IbP
    self.A_vf.add(( -1j*iDot(nbc_tens,fluc_rhou*iConj(X)) ).ufl_tens*coord.J_hat*self.all_ds)

    # ------------------------ BC term for Input/Output analysis
    if param.Case.AnalysisMode in ['Input-Output']:
        # Iterate through all boundaries, at which forcing is applied
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            # First subtract the boundary term from advection
            # Correction wrt to index notation: We need to remove the rho*u term, not just the u! 
            self.A_vf.add(( 1j*iDot(nbc_tens,fluc_rhou*iConj(X)) ).ufl_tens*coord.J_hat*self.ds(boundary_index))
            # Then add the forcing at the boundary
            self.A_vf.add(( -1*iDot(nbc_tens,forcing_tens)*iConj(X) ).ufl_tens*coord.J_hat*self.ds(boundary_index))
                
