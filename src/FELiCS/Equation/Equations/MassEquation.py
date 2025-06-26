from ufl import dx
from FELiCS.Misc.tensorUtils import (
    iGrad,
    iDot,
    iConj
)

from .EquationTemplate      import EquationTemplate
from FELiCS.Misc.logging    import Logger

# Get the logger
logger = Logger.get_logger("felics")

class MassEquation(EquationTemplate):
    """
    Class representing the mass conservation equation.

    Parameters
    ----------
    eqColl : EquationCollection
        The equation collection object.
    fluc : Fluctuations
        The fluctuations object.
    X : Function
        The solution function.
    param : Parameters
        The parameters object.

    Notes
    -----
    This class inherits from EquationTemplate and implements the specific
    methods for the mass conservation equation.
    """

    def __init__(self, eqColl, fluc, X, param):
        """
        Initialize the MassEquation object.

        Parameters
        ----------
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuations object.
        X : Function
            The solution function.
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

        Parameters
        ----------
        weakForm : Form
            The weak form object.
        mean : Function
            The mean function.

        Notes
        -----
        This function adds the time derivative term to the weak form,
        only if density fluctuations are considered.
        """
        # ------------------------ Time derivative term used
        # Only if density fluctuations are considered
        if 'rho' in self.param.getTransportedQuantityList():
            weakForm += (self.fluc.rho * iConj(self.X)).ufl_tens * self.J_hat * dx

    def addLinearExpression(self, weakForm, mean):
        """
        Add the linear expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object.
        mean : Function
            The mean function.

        Notes
        -----
        This function builds the weak form of the linearized mass conservation equation,
        in tensorial framework.
        """
        # ------------------------ Advection terms
        # The advection term is integrated by parts
        # Volume term from IbP
        weakForm.add((1j * iDot(iGrad(iConj(self.X)), self.fluc.rhou)).ufl_tens * self.J_hat * dx)
        # Boundary term from IbP
        weakForm.add((-1j * iDot(self.n, self.fluc.rhou * iConj(self.X))).ufl_tens * self.J_hat * self.all_ds)

        # ------------------------ BC term for Input/Output analysis
        if self.param.Case.AnalysisMode in ['Input-Output']:
            # Iterate through all boundaries, at which forcing is applied
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:
                # First subtract the boundary term from advection
                # Correction wrt to index notation: We need to remove the rho*u term, not just the u!
                weakForm.add((1j * iDot(self.n, self.fluc.rhou * iConj(self.X))).ufl_tens * self.J_hat * self.ds(boundary_index))
                # Then add the forcing at the boundary
                weakForm.add((-1 * iDot(self.n, mean.u_forcing * mean.rho) * iConj(self.X)).ufl_tens * self.J_hat * self.ds(boundary_index))


    def addNonlinearExpression(self, weakForm, mean):
        """
        Add the nonlinear expression to the weak form.

        Notes
        -----
        This function builds the weak form of the nonlinear
        mass conservation equation, in tensorial framework.
        """

        # ------------------------ Advection terms
        # The advection term is integrated by parts
        # Volume term from IbP
        weakForm.add((  1j * iDot(iGrad(iConj(self.X)),mean.rho*mean.u)).ufl_tens * self.J_hat * dx)
        # Boundary term from IbP
        weakForm.add(( -1j * iDot(self.n,mean.rho*mean.u * iConj(self.X)) ).ufl_tens * self.J_hat * self.all_ds)
    



