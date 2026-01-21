#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
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

    The sponge term applies damping to fluctuations or deviations from target
    profiles in computational domains. It stabilizes numerical simulations by
    attenuating reflections and enforcing desired flow behavior in designated
    regions.

    **Initialize the SpongeTerm object**

    Parameters
    ----------
    eqColl : EquationCollection
        The equation collection object.
    fluc : Fluctuations
        The fluctuations object.
    X : list
        The list of solution variables.
    param : Parameters
        The parameters object.

    Attributes
    ----------
    fluc : Fluctuations
        Fluctuating fields used in the formulation.
    X : list of Function
        Trial/test functions for each variable.
    param : Parameters
        Configuration and problem parameters.
    J_hat : Expression
        Jacobian determinant for integration.
    all_ds : Measure
        Boundary integration measure.
    n : FacetNormal
        Unit normal vector on boundaries.

    Notes
    -----
    Discontinuous Galerkin schemes are not supported in this tensorial framework.
    """

    def __init__(self, index, eqColl, fluc, X, param):
        """
        Initialize the SpongeTerm object.

        Parameters
        ----------
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
        super().__init__(index, eqColl, fluc, X, param)


    def addWeightMatrixExpression(self, weakForm, mean):
        """
        Add the weight matrix expression to the weak form.

        Parameters
        ----------
        weakForm : ufl.Form
            The weak form object.
        mean : MeanField
            The mean field object.

        Notes
        -----
        This method does not apply any sponge term contributions.
        """
        # nothing to add for the sponge term
        pass

    def addLinearExpression(self, weakForm, mean):
        """
        Add the linear expression to the weak form.

        Parameters
        ----------
        weakForm : ufl.Form
            The weak form object.
        mean : MeanField
            The mean field object.

        Notes
        -----
        Applies linear sponge damping terms to the fluctuation variables,
        as specified in the parameter set. For velocity variables, a dot
        product is applied. Each variable's sponge term is weighted by the
        sponge strength field `spg`.
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

        Parameters
        ----------
        weakForm : ufl.Form
            The weak form object.
        mean : MeanField
            The mean field object.

        Notes
        -----
        Applies nonlinear sponge damping terms based on the deviation of
        mean flow variables from their designated target profiles. These
        target values are retrieved dynamically and matched to each variable.
        Velocity terms are handled with tensor dot products.
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
                    target_u = mean._fieldDict['u_target'].getTensor()
                    weakForm.add(( -1j*mean.spg*iDot(mean_var-target_u,iConj(X[varNum])) ).ufl_tens*J_hat*dx)
                else: 
                    # for every other name there has to exist a field in the mean flow dictionary with the name and the suffix '_target'
                    # TODO: what to do if the field does not exist? Logging: throw error
                    target_name = varID+"_target"
                    target_tens = mean._fieldDict[target_name].getTensor()
                    weakForm.add(( -1j*mean.spg*(mean_var-target_tens)*iConj(X[varNum]) ).ufl_tens*J_hat*dx)
