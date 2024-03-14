from ufl import dx
from FELiCS.tensorUtils import (
    iDot,
    iConj
)
from FELiCS.functions import printError, printDebug

def addSpongeEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the sponge
    term in tensorial framework.
    '''

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')        

    J_hat = self._coordinateSystem.J_hat
    
    # Looping over all linearized equations
    for eqID in param.Case.SetOfEquations :
        if not (param.Case.SetOfEquations[eqID]['Equation'] == 'None' or \
            param.Case.SetOfEquations[eqID]['Variable'] == 'None'):
            
            #  Get the corresponding variable
            varID = param.Case.SetOfEquations[eqID]['Variable']
            varNum = param.SolutionList.index(varID)
            printDebug(True, "-- -> Sponge term for %s-fluc -> X[%d]." % (varID,varNum))
            
            # Dynamically get the corresponding fluctuation field
            fluc_var = getattr(fluc, '%s' % varID)
            
            # Apply the sponge
            if varID == 'u': # For u we need the dot product with X
                self.A_vf.add(( -1j*mean.spg*iDot(fluc_var,iConj(X[varNum])) ).ufl_tens*J_hat*dx)
            else:
                self.A_vf.add(( -1j*mean.spg*fluc_var*iConj(X[varNum]) ).ufl_tens*J_hat*dx)
        

    # # Assuming velocity fluctuations are ALWAYS considered
    # self.A_vf.add(( -1j*mean.spg*iDot(fluc.u,iConj(X[0])) ).ufl_tens*J_hat*dx)
    
    # # Looping over all other linearized equation    
    # id_p = param.SolutionList.index('p')
    # self.A_vf.add(( -1j*mean.spg*fluc.p*iConj(X[id_p]) ).ufl_tens*J_hat*dx)
    # # ------------------------ Compressible
    # # --> Assuming rho is the extra state variable!
    # if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
    #     printDebug(True, '-- -> Adding sponge damping for density fluctuations.')
    #     id_rho = param.SolutionList.index('rho')
    #     self.A_vf.add(( -1j*mean.spg*fluc.rho*iConj(X[id_rho]) ).ufl_tens*J_hat*dx)
