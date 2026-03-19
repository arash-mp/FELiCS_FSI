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
# Standard libraries

# Third party libraries
from ufl import dx

# Local Libraries and methods
from FELiCS.Equation.Boundary                   import BoundaryType
from FELiCS.Equation.Equations.EquationTemplate import EquationTemplate
from FELiCS.Misc.logging                        import Logger, log_and_raise
from FELiCS.Misc.tensorUtils                    import (
    i_conj,
    i_div,
    i_dot,
    i_grad,
)

# Get the logger
logger = Logger.get_logger("felics")


class SpeciesEquation(EquationTemplate):
    """
    Class representing the species transport equation.

    This class formulates the species transport equation in a tensorial framework,
    capturing the effects of advection, diffusion, and chemical reactions. It supports
    boundary forcing and input-output analysis, while enforcing integration by parts
    to enable more robust numerical handling of fluxes and boundary conditions.

    **Initialize the SpeciesEquation object**

    Parameters
    ----------
    eqColl : EquationCollection
        The equation collection object.
    fluc : Fluctuations
        The fluctuations object.
    X : Function
        The trial/test function in the weak formulation.
    species : str
        The species name being transported.
    param : Parameters
        The configuration and simulation parameters.

    Attributes
    ----------
    species : str
        Name of the species for which the equation is formulated.
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
    Discontinuous Galerkin schemes are not supported in this tensorial framework.
    """


    def __init__(
        self,
        index,
        eqColl,
        fluc,
        X,
        species,
        param
    ):
        """
        Initialize the SpeciesEquation class.

        Parameters
        ----------
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuation object.
        X : Function
            The function representing the mesh coordinates.
        species : str
            The species name.
        param : Parameters
            The parameters object.
        """

        # Disclaimer
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

        self.species = species


    def add_weight_matrix_expression(
        self,
        weakForm,
        mean
    ):
        """
        Add the weight matrix expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object.
        mean : MeanFlow
            The mean flow object.

        Notes
        -----
        Adds the time derivative term to the weak form for the species equation.
        This term incorporates fluctuations weighted by the mean density.
        """

        # Time derivative term
        weakForm.add((self.fluc.y(self.species) * i_conj(self.X) * mean.rho).ufl_tens * self.J_hat * dx)

    def add_nonlinear_expression(self):
        """
        Add the nonlinear expression to the weak form.

        Warning
        -----
        This method is currently a placeholder and not implemented.
        """

        pass

    def add_linear_expression(
        self,
        weakForm,
        mean
    ):
        """
        Construct the weak form of the linearized species transport equation.

        Parameters
        ----------
        weakForm : Form
            The weak form object.
        mean : MeanFlow
            The mean flow object.

        Notes
        -----
        Constructs the weak form using a tensorial formulation. The method
        includes:
        - Integration by parts for advection terms to capture boundary contributions.
        - Volume-only integration of diffusion terms, effectively imposing Neumann boundary conditions.
        - Reaction terms based on a KaiserCnF2023 mechanism if specified.
        - Forcing terms for input-output analysis, including both body and boundary forcing.
        A warning is issued if the case parameter `m > 0`, as it has not been validated.
        """
       
        param   = self.param
        J_hat   = self.J_hat
        species = self.species
        X       = self.X
        fluc    = self.fluc

        if not param.Case.m == 0:
            logger.warning('Species eq. with m > 0 not validated yet for tensor. Treat results with care.')
        
        # ----------------------------------------- Advection term
        # This term is integrated by parts
        # Volume term from IbP
        # The volume term seems to introduce a small error (~1e-12) in cartesian nates wrt. previous implementation
        ibp = True
        if ibp:
            logger.debug("Using integration by parts for advection term.")
            weakForm.add(( 1j*fluc.y(species)*i_div(mean.rho*mean.u*i_conj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( 1j*mean.y(species)*i_div(mean.rho*fluc.u*i_conj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( 1j*mean.y(species)*i_div(fluc.rho*mean.u*i_conj(X)) ).ufl_tens*J_hat*dx)

            # Boundary term from IbP
            weakForm.add(
                ( -1j * i_dot(
                    self.n,
                    fluc.y(species) * mean.rho * mean.u * i_conj(X),
                )).ufl_tens * J_hat * self.all_ds
            )

            weakForm.add(
                ( -1j * i_dot(
                    self.n,
                    mean.y(species) * mean.rho * fluc.u * i_conj(X),
                )).ufl_tens * J_hat * self.all_ds
            )

            weakForm.add(
                ( -1j * i_dot(
                    self.n,
                    mean.y(species) * fluc.rho * mean.u * i_conj(X),
                )).ufl_tens * J_hat * self.all_ds
            )

        else:
            logger.debug("NOT using integration by parts for advection term.")

            weakForm.add(
                ( -1j * i_dot(
                i_grad(fluc.y(species)),
                mean.rho * mean.u * i_conj(X),
                )).ufl_tens * J_hat * dx
            )

            weakForm.add(
                ( -1j * i_dot(
                    i_grad(mean.y(species)),
                    mean.rho * fluc.u * i_conj(X),
                )).ufl_tens * J_hat * dx
            )

            weakForm.add(
                ( -1j * i_dot(
                    i_grad(mean.y(species)),
                    fluc.rho * mean.u * i_conj(X),
                )).ufl_tens * J_hat * dx
            )
            
    
    
        # ----------------------------------------- Diffusion term
        # NOTE: similarly to what is done in the mom. eq., the diffusion term is integrated by parts but only the 
        # volume part is added to the eqs. --> neglecting the boundary part allows to set a Neumann condition  
        weakForm.add(
            ( -1j * mean.d(species) * (i_dot(
                i_grad(fluc.y(species)),
                i_grad(i_conj(X)),
            ))).ufl_tens * J_hat * dx
        )

        weakForm.add(
            ( -1j * fluc.d(species) * (i_dot(
                i_grad(mean.y(species)),
                i_grad(i_conj(X)),
            ))).ufl_tens * J_hat * dx
        )
        
        #--------------------------------Reaction
        if self.param.Mixture.reaction_mechanism['type'] == 'KaiserCnF2023':
            reaction = mean.rr_prefactor * mean.rho * (fluc.y(species) - 2 * fluc.y(species) * mean.y(species))\
                            + mean.rr_prefactor * fluc.rho * (mean.y(species) - mean.y(species) * mean.y(species))
            weakForm.add((1j * reaction * i_conj(X)).ufl_tens * J_hat * dx)


        # ----------------------------------------- Input/Output forcing
        if self.param.Case.AnalysisMode in ['Input-Output']:

            # add body forcing
            if param.IOResolvent.ForcingMode == 'Body':
                weakForm.add(
                    ( mean.forcing(species) * i_conj(X)).ufl_tens * J_hat * dx
                )

            # Iterate through all boundaries, at which boundary forcing is applied
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:

                # First subtract the part added in a few lines above...
                weakForm.add(
                    ( 1j * i_dot(
                        self.n,
                        fluc.y(species) * mean.rho * mean.u * i_conj(X),
                    )).ufl_tens * J_hat * self.ds(boundary_index)
                )

                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                weakForm.add(
                    ( -1j * i_dot(
                        self.n,
                        mean.u * mean.forcing(species) * i_conj(X),
                    ) ).ufl_tens * J_hat * self.ds(boundary_index)
                )
                    
