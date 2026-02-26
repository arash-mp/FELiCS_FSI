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
# Third party libraries
from ufl import dx


# Local Libraries and methods
from FELiCS.Equation.Equations.EquationTemplate import EquationTemplate
from FELiCS.Misc.logging                        import Logger, log_and_raise
from FELiCS.Misc.tensorUtils                    import (
    i_conj,
    i_dot,
    i_grad,
)


# Get the logger
logger = Logger.get_logger("felics")

class MassEquation(EquationTemplate):
    """
    Class representing the mass conservation equation.

    This class implements the weak form expressions for the mass conservation
    equation, both linear and nonlinear, in the tensorial framework. It
    handles the addition of relevant terms such as advection and boundary
    contributions, and accounts for specific analysis modes like
    Input-Output.

    **Initialize the MassEquation object**

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

    Attributes
    ----------
    fluc : Fluctuations
        Fluctuating quantities used in the formulation.
    X : Function
        Trial/test function in the variational formulation.
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
    This class inherits from EquationTemplate and raises an exception if
    a Discontinuous Galerkin scheme is selected, which is not implemented
    in the tensorial framework.
    """


    def __init__(
        self,
        index,
        eqColl,
        fluc,
        X,
        param,
    ):
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
            log_and_raise(logger, 'Discontinuous Galerkin not implemented in tensorial framework.', Exception)
        # initialize variables in template class
        super().__init__(
            index,
            eqColl,
            fluc,
            X,
            param,
        )

    def add_weight_matrix_expression(
        self,
        weakForm,
        mean,
    ):
        """
        Add the weight matrix expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the expression is added.
        mean : Function
            The mean function providing averaged quantities.

        Notes
        -----
        Adds the time-derivative term to the weak form if density is among the
        transported quantities.
        """
        # ------------------------ Time derivative term used
        # Only if density fluctuations are considered
        if 'rho' in self.param.get_transported_quantity_list():
            weakForm += (self.fluc.rho * i_conj(self.X)).ufl_tens * self.J_hat * dx

    def add_linear_expression(
        self,
        weakForm,
        mean,
    ):
        """
        Add the linear expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the expression is added.
        mean : Function
            The mean function providing averaged quantities.

        Notes
        -----
        Constructs the weak form of the linearized mass conservation equation,
        including volume and boundary terms via integration by parts. Special
        handling is added for boundary forcing in Input-Output analysis mode.
        """
        # ------------------------ Advection terms
        # The advection term is integrated by parts

        # Volume term from IbP
        weakForm.add((1j * i_dot(
                i_grad(i_conj(self.X)),
                self.fluc.rhou,
            )).ufl_tens * self.J_hat * dx
        )

        # Boundary term from IbP
        weakForm.add(
            (-1j * i_dot(
                self.n, 
                self.fluc.rhou * i_conj(self.X),
            )).ufl_tens * self.J_hat * self.all_ds
        )

        # ------------------------ BC term for Input/Output analysis
        if self.param.Case.AnalysisMode in ['Input-Output']:
            # Iterate through all boundaries, at which forcing is applied
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:
                # First subtract the boundary term from advection
                weakForm.add((1j * i_dot(
                    self.n,
                    self.fluc.u * mean.rho * i_conj(self.X),
                )).ufl_tens * self.J_hat * self.ds(boundary_index))

                # Then add the forcing at the boundary
                weakForm.add(
                    (-1 * i_dot(
                        self.n,
                        mean.u_forcing * mean.rho,
                    ) * i_conj(self.X)).ufl_tens * self.J_hat * self.ds(boundary_index)
                )


    def add_nonlinear_expression(
        self,
        weakForm,
        mean,
    ):
        """
        Add the nonlinear expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the expression is added.
        mean : Function
            The mean function providing averaged quantities.

        Notes
        -----
        Constructs the weak form of the full nonlinear mass conservation
        equation using integration by parts for advection terms.
        """

        # ------------------------ Advection terms
        # The advection term is integrated by parts
        # Volume term from IbP
        weakForm.add(
            (1j * i_dot(
                i_grad(i_conj(self.X)),
                mean.rho*mean.u,
            )).ufl_tens * self.J_hat * dx
        )

        # Boundary term from IbP
        weakForm.add(
            (-1j * i_dot(
                self.n,
                mean.rho*mean.u * i_conj(self.X),
            )).ufl_tens * self.J_hat * self.all_ds
        )
    