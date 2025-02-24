from    ufl                     import dx
from    .EquationTemplate       import EquationTemplate
from    FELiCS.Misc.logging     import Logger
from    FELiCS.Misc.tensorUtils import (
    Tensor,
    iDot,
    iConj
)

# Get the logger
logger = Logger.get_logger("felics")

class SpongeTerm(EquationTemplate):
    """
    Class representing the sponge term in the equation.
    """

    def __init__(self, eqColl, fluc, X, param):
        """
        Initialize the SpongeTerm object.

        Parameters:
        -----------
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuations object.
        X : list
            The list of solution variables.
        param : Parameters
            The parameters object.
        """
        # Disclaimers
        if param.Numerics.NumericalScheme in ['Discontinuous Galerkin']:
            logger.error('Discontinuous Galerkin not implemented in tensorial framework.')
            raise Exception('Discontinuous Galerkin not implemented in tensorial framework.')        

        # initialize variables in template class
        super().__init__(eqColl, fluc, X, param)


    def addWeightMatrixExpression(self, weakForm, mean):
        """
        Add the weight matrix expression to the weak form.

        Parameters:
        -----------
        weakForm : ufl.Form
            The weak form object.
        mean : MeanField
            The mean field object.
        """
        # nothing to add
        pass

    def addLinearExpression(self, weakForm, mean):
        """
        Add the linear expression to the weak form.

        Parameters:
        -----------
        weakForm : ufl.Form
            The weak form object.
        mean : MeanField
            The mean field object.
        """
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
                varNum = param.Case.SolutionList.index(varID)
                logger.debug("Adding sponge term for %s-fluc: X[%d]." % (varID,varNum))
                
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
        # id_p = param.Case.SolutionList.index('p')
        # weakForm.add(( -1j*mean.spg*fluc.p*iConj(X[id_p]) ).ufl_tens*J_hat*dx)
        # # ------------------------ Compressible
        # # --> Assuming rho is the extra state variable!
        # if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
        #     printDebug(True, '-- -> Adding sponge damping for density fluctuations.')
        #     id_rho = param.Case.SolutionList.index('rho')
        #     weakForm.add(( -1j*mean.spg*fluc.rho*iConj(X[id_rho]) ).ufl_tens*J_hat*dx)

    def addNonlinearExpression(self, weakForm, mean):
        """
        Add the nonlinear expression to the weak form.
        """
 
        J_hat = self.J_hat
        X     = self.X
        param = self.param
        
        # Looping over all linearized equations
        for eqID in param.Case.SetOfEquations :
            if not (param.Case.SetOfEquations[eqID]['Equation'] == 'None' or \
                param.Case.SetOfEquations[eqID]['Variable'] == 'None'):
                
                #  Get the corresponding variable
                varID = param.Case.SetOfEquations[eqID]['Variable']
                varNum = param.SolutionList.index(varID)
                logger.debug("Adding sponge term for %s-fluc: X[%d]." % (varID,varNum))
                
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
            


