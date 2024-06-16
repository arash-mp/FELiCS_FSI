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

class SpongeTerm(EquationTemplate):

    def __init__(self,eqColl,fluc,X,param):
 
        # Disclaimers
        if param.NumericalScheme in ['Discontinuous Galerkin']:
            printError('Discontinuous Galerkin not implemented in tensorial framework.')        

        # initialize variables in template class
        super().__init__(eqColl,fluc,X,param)


    def addWeightMatrixExpression(self, weakForm, mean):
        # nothing to add
        pass

    def addLinearExpression(self, weakForm, mean):
        '''
        This function builds the weak form of the sponge
        term in tensorial framework.
        '''
    
        J_hat = self.J_hat
        fluc  = self.fluc
        X     = self.X
        param = self.param
        
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
                    weakForm.add(( -1j*mean.spg*iDot(fluc_var,iConj(X[varNum])) ).ufl_tens*J_hat*dx)
                else:
                    weakForm.add(( -1j*mean.spg*fluc_var*iConj(X[varNum]) ).ufl_tens*J_hat*dx)
            
    
        # # Assuming velocity fluctuations are ALWAYS considered
        # weakForm.add(( -1j*mean.spg*iDot(fluc.u,iConj(X[0])) ).ufl_tens*J_hat*dx)
        
        # # Looping over all other linearized equation    
        # id_p = param.SolutionList.index('p')
        # weakForm.add(( -1j*mean.spg*fluc.p*iConj(X[id_p]) ).ufl_tens*J_hat*dx)
        # # ------------------------ Compressible
        # # --> Assuming rho is the extra state variable!
        # if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
        #     printDebug(True, '-- -> Adding sponge damping for density fluctuations.')
        #     id_rho = param.SolutionList.index('rho')
        #     weakForm.add(( -1j*mean.spg*fluc.rho*iConj(X[id_rho]) ).ufl_tens*J_hat*dx)

    def addNonlinearExpression(self, weakForm, mean):
 
        J_hat = self.J_hat
        fluc  = self.fluc
        X     = self.X
        param = self.param
        
        # Looping over all linearized equations
        for eqID in param.Case.SetOfEquations :
            if not (param.Case.SetOfEquations[eqID]['Equation'] == 'None' or \
                param.Case.SetOfEquations[eqID]['Variable'] == 'None'):
                
                #  Get the corresponding variable
                varID = param.Case.SetOfEquations[eqID]['Variable']
                varNum = param.SolutionList.index(varID)
                printDebug(True, "-- -> Sponge term for %s-fluc -> X[%d]." % (varID,varNum))
                
                # Dynamically get the corresponding fluctuation field
                mean_var = getattr(mean, '%s' % varID)

                # Apply the sponge
                if varID == 'u': # For u we need the dot product with X
                    target_u = Tensor(mean._fieldDict['u_target'], self.coordinateSystem)
                    weakForm.add(( -1j*mean.spg*iDot(mean_var-target_u,iConj(X[varNum])) ).ufl_tens*J_hat*dx)
                elif varID == 'p':
                    target_p = Tensor(mean._fieldDict['p_target'], self.coordinateSystem)
                    weakForm.add(( -1j*mean.spg*(mean_var-target_p)*iConj(X[varNum]) ).ufl_tens*J_hat*dx)
                else:
                    weakForm.add(( -1j*mean.spg*mean_var*iConj(X[varNum]) ).ufl_tens*J_hat*dx)
            


