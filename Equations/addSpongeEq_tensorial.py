from ufl import dx
from tensorUtils import (
    iDot,
    iConj
)
from functions import printError

def addSpongeEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the sponge
    term in tensorial framework.
    '''

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')        

    J_hat = self._coordinateSystem.J_hat

    # ------------------------ Incompressible
    self.A_vf.add(( -1*mean.spg*iDot(fluc.u,iConj(X[0])) ).ufl_tens * J_hat * dx)
    id_p = param.SolutionList.index('p')
    self.A_vf.add(( -1*mean.spg*fluc.p*iConj(X[id_p]) ).ufl_tens * J_hat * dx)
    # self.A_vf.add(( -1*mean.spg * (iDot(fluc.u,iConj(X[0])) + fluc.p*iConj(X[id_p])) ).ufl_tens * J_hat * dx)
    
    # ------------------------ Compressible
    # --> Assuming rho is the extra state variable!
    id_rho = param.SolutionList.index('rho')
    if not param.Case.SetOfEquations['Energy']['Equation'] == 'none':
        self.A_vf.add(( -1*mean.spg * fluc.rho*iConj(X[id_rho]) ).ufl_tens * J_hat * dx)
