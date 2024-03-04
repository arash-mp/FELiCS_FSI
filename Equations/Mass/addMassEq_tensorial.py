from ufl import dx
from tensorUtils import (
    Tensor,
    as_vector,
    iGrad,
    iDot,
    iConj
)
from functions import printDebug, printError

def addMassEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    mass conservation equation, in tensorial framework.
    '''

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')        

    J_hat = self._coordinateSystem.J_hat

    # ------------------------ Time derivative term used
    # Only if density fluctuations are considered
    if 'rho' in param.Case.getTransportedQuantityList():
        self.B_vf.add( (fluc.rho * iConj(X)).ufl_tens * J_hat * dx)

    # ------------------------ Advection terms
    # The advection term is integrated by parts
    # Volume term from IbP
    self.A_vf.add((  1j * iDot(iGrad(iConj(X)),fluc.rhou)).ufl_tens * J_hat * dx)
    # Boundary term from IbP
    self.A_vf.add(( -1j * iDot(self.n,fluc.rhou * iConj(X)) ).ufl_tens * J_hat * self.all_ds)

    # ------------------------ BC term for Input/Output analysis
    if param.Case.AnalysisMode in ['Input-Output']:
        # Iterate through all boundaries, at which forcing is applied
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            # First subtract the boundary term from advection
            # Correction wrt to index notation: We need to remove the rho*u term, not just the u! 
            self.A_vf.add(( 1j * iDot(self.n,fluc.rhou * iConj(X)) ).ufl_tens * J_hat * self.ds(boundary_index))
            # Then add the forcing at the boundary
            self.A_vf.add(( -1 * iDot(self.n,mean.u_forcing* mean.rho) * iConj(X) ).ufl_tens * J_hat * self.ds(boundary_index))
