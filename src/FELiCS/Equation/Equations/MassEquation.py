from ufl import dx
from FELiCS.Misc.tensorUtils import (
    Tensor,
    as_vector,
    iGrad,
    iDot,
    iConj
)
from FELiCS.Misc.functions import printDebug, printError

from .EquationTemplate import EquationTemplate

class MassEquation(EquationTemplate):

    def __init__(self,eqColl,fluc,X,param):
 
        # Disclaimers
        if param.NumericalScheme in ['Discontinuous Galerkin']:
            printError('Discontinuous Galerkin not implemented in tensorial framework.')        

        # initialize variables in template class
        super().__init__(eqColl,fluc,X,param)


    def addWeightMatrixExpression(self, weakForm, mean):
 
        # ------------------------ Time derivative term used
        # Only if density fluctuations are considered
        if 'rho' in self.param.Case.getTransportedQuantityList():
            weakForm.add( (self.fluc.rho * iConj(self.X)).ufl_tens * self.J_hat * dx)
    

    def addLinearExpression(self, weakForm, mean):

        '''
        This function builds the weak form of the linearized
        mass conservation equation, in tensorial framework.
        '''
              # ------------------------ Advection terms
        # The advection term is integrated by parts
        # Volume term from IbP
        weakForm.add((  1j * iDot(iGrad(iConj(self.X)),self.fluc.rhou)).ufl_tens * self.J_hat * dx)
        # Boundary term from IbP
        weakForm.add(( -1j * iDot(self.n,self.fluc.rhou * iConj(self.X)) ).ufl_tens * self.J_hat * self.all_ds)
    
        # ------------------------ BC term for Input/Output analysis
        if self.param.Case.AnalysisMode in ['Input-Output']:
            # Iterate through all boundaries, at which forcing is applied
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:
                # First subtract the boundary term from advection
                # Correction wrt to index notation: We need to remove the rho*u term, not just the u! 
                weakForm.add(( 1j * iDot(self.n,self.fluc.rhou * iConj(self.X)) ).ufl_tens * self.J_hat * self.ds(boundary_index))
                # Then add the forcing at the boundary
                weakForm.add(( -1 * iDot(self.n,mean.u_forcing* mean.rho) * iConj(self.X) ).ufl_tens * self.J_hat * self.ds(boundary_index))


    def addNonlinearExpression(self):
        pass


