from ufl import dx
from tensorUtils import (
    iDot,
    iConj
)
from functions import printError, printDebug

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
    self.A_vf.add(( -1j*mean.spg*iDot(fluc.u,iConj(X[0])) ).ufl_tens*J_hat*dx)
    id_p = param.SolutionList.index('p')
    self.A_vf.add(( -1j*mean.spg*fluc.p*iConj(X[id_p]) ).ufl_tens*J_hat*dx)
    # ------------------------ Compressible
    # --> Assuming rho is the extra state variable!
    if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
        printDebug(True, '-- Adding sponge damping for density fluctuations.')
        id_rho = param.SolutionList.index('rho')
        self.A_vf.add(( -1j*mean.spg*fluc.rho*iConj(X[id_rho]) ).ufl_tens*J_hat*dx)
