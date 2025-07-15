from ufl import (
    dx
)
from FELiCS.Misc.tensorUtils import (
    iDot,
    iDiv,
    iGrad,
    iConj
)
from    .EquationTemplate   import EquationTemplate
from    FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")


class EnthalpyEquation(EquationTemplate):
    """
    Class representing the enthalpy conservation equation.

    This class formulates the enthalpy conservation equation using a 
    tensorial framework. It extends the `EquationTemplate` and implements 
    methods for adding weight matrices, linear expressions, and nonlinear 
    expressions to the weak form.

    **Initialize the EnthalpyEquation object**

    Parameters
    ----------
    index : int
        The index of the equation in the equation collection.
    eqColl : EquationCollection
        The equation collection object.
    fluc : Fluctuations
        The fluctuations object containing fluctuating fields.
    X : Function
        The function representing the mesh coordinates.
    param : Parameters
        The parameters object used to configure the numerical scheme.

    Attributes
    ----------
    J_hat : float or Function
        Jacobian determinant or related geometric scaling factor.
    X : Function
        The mesh coordinate function.
    fluc : Fluctuations
        Container of fluctuation fields.
    all_ds : Measure
        Surface measure used for boundary integrals.

    Notes
    -----
    If the numerical scheme specified in the parameters is 
    'Discontinuous Galerkin', an error is raised because that scheme 
    is not supported in the tensorial framework.
    """

    def __init__(self, index, eqColl, fluc, X, param):
        """
        Initialize the EnthalpyEquation object.

        Parameters
        ----------
        index : int
            The index of the equation in the equation collection.
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuations object.
        X : Function
            The function representing the mesh coordinates.
        param : Parameters
            The parameters object.
        
        Notes
        -----
        - If the numerical scheme is 'Discontinuous Galerkin', an error is raised
          because it is not implemented in the tensorial framework.
        """
        # Disclaimer
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
        weakForm : Form
            The weak form object where the equation terms are added.
        mean : MeanFields
            The mean fields object containing averaged field variables.

        Notes
        -----
        This method adds the time derivative terms to the weak form using the
        tensorial representation of fluctuating enthalpy and pressure.
        """
        #  Time derivative terms
        weakForm.add((mean.rho * self.fluc.h * iConj(self.X)).ufl_tens * self.J_hat * dx)
        weakForm.add((-1 * self.fluc.p * iConj(self.X)).ufl_tens * self.J_hat * dx)
    
     
    def addNonlinearExpression(self):
        """
        Add the nonlinear expression to the weak form.

        Notes
        -----
        This method is a placeholder and currently does not contribute any 
        nonlinear terms to the weak form.
        """
        pass

    def addLinearExpression(self, weakForm, mean):
        """
        Construct the weak form of the linearized enthalpy conservation equation.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the linearized terms are added.
        mean : MeanFields
            The mean fields object containing averaged field variables.

        Notes
        -----
        Constructs the linearized weak form by incorporating:
        - Volume integrals of advection terms remaining after partial integration.
        - Boundary integrals resulting from the partial integration of advection terms.
        - Diffusion terms using gradient and dot product operations.

        Debugging messages may be emitted if enabled via the parameters.
        """

        J_hat   = self.J_hat
        X       = self.X
        fluc    = self.fluc
        
            
        # ------------------------  Advection terms
        # Add volume integral of advection terms that remain after partial integration:
        weakForm.add((1j * iDiv(iConj(X) * mean.rho * mean.u) * fluc.h).ufl_tens * J_hat * dx)
        weakForm.add((1j * iDiv(iConj(X) * fluc.rho * mean.u) * mean.he).ufl_tens * J_hat * dx)
        weakForm.add((1j * iDiv(iConj(X) * mean.rho * fluc.u) * mean.he).ufl_tens * J_hat * dx)
        # Add boundary integrals resulting from said partial integration:
        weakForm.add((-1j * iDot(mean.rho * mean.u * fluc.h * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
        weakForm.add((-1j * iDot(fluc.rho * mean.u * mean.he * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
        weakForm.add((-1j * iDot(mean.rho * fluc.u * mean.he * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
            
        # ------------------------  Diffusion terms
        weakForm.add((-1j * iDot(mean.alpha * iGrad(fluc.h), iGrad(iConj(X)))).ufl_tens * J_hat * dx)
        weakForm.add((-1j * iDot(fluc.alpha * iGrad(mean.he), iGrad(iConj(X)))).ufl_tens * J_hat * dx)
