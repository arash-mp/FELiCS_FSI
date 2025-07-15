from ufl import (
    dx
)
from FELiCS.Misc.tensorUtils import (
    Tensor,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
)
from .EquationTemplate      import EquationTemplate
from FELiCS.Misc.logging    import Logger

# Get the logger
logger = Logger.get_logger("felics")

class MomentumEquation(EquationTemplate):
    """
    Class representing the momentum conservation equation.

    This class formulates the momentum conservation equation in a tensorial framework.
    It includes convective, pressure gradient, and diffusion terms while supporting 
    both weak and strong formulations. The class integrates with the overall equation 
    collection and handles interactions with mean fields and fluctuations.

    **Initialize the MomentumEquation object**

    Parameters
    ----------
    index : int
        Index of the equation in the system.
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
    If the numerical scheme is 'Discontinuous Galerkin', an error will be raised 
    since it is not implemented in the tensorial framework.
    """

    def __init__(self, index, eqColl, fluc, X, param):
        """
        Initialize the MomentumEquation instance.

        Parameters
        ----------
        index : int
            Index of the equation in the system.
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
        If the numerical scheme is 'Discontinuous Galerkin', an error will be raised since it is not implemented in the tensorial framework.
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
            The weak form object.
        mean : MeanFields
            The mean fields object.

        Notes
        -----
        This method incorporates time derivative terms into the weak form for the momentum equation.
        """

        # Time derivative term
        weakForm += (iDot(mean.rho*self.fluc.u, iConj(self.X))).ufl_tens*self.J_hat*dx

    def addLinearExpression(self,weakForm,mean):
        """
        Construct the weak form of the linearized momentum conservation equation.

        Parameters
        ----------
        weakForm : Form
            The weak form object.
        mean : MeanFields
            The mean fields object.

        Notes
        -----
        - This function builds the weak form of the momentum conservation equation 
          in a tensorial framework.
        - Convective, pressure gradient, and diffusion terms are considered.
        - Integration by parts is optionally applied depending on the coordinate system.
        - Certain terms may require adjustments in future implementations.
        """
      
        fluc  = self.fluc
        X     = self.X
        J_hat = self.J_hat

        # ------------------------ Convective terms
        int_by_parts = False
        if int_by_parts and self.param.Case.CoordinateSystem=='Cartesian':
            logger.debug(" -> Using integration by parts for convection term.")
            # Volume term from integration by parts
            weakForm += ( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*mean.u)),fluc.u) ).ufl_tens*J_hat*dx
            weakForm += ( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*fluc.u)),mean.u) ).ufl_tens*J_hat*dx
            weakForm += ( 1j*iDot(iDiv(iOuter(iConj(X),fluc.rho*mean.u)),mean.u) ).ufl_tens*J_hat*dx
            # Boundary term from integration by parts
            weakForm += ( -1j*mean.rho*iDot(iDot(iOuter(fluc.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds
            weakForm += ( -1j*mean.rho*iDot(iDot(iOuter(mean.u,iConj(X)),fluc.u),self.n) ).ufl_tens*J_hat*self.all_ds
            weakForm += ( -1j*fluc.rho*iDot(iDot(iOuter(mean.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds
            #elif self.param.Case.CoordinateSystem =='Cylindrical':
            #    # In cyl , a singular term error arise for the boundary term in the tensor framework
            #    # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
            #    # This should be fixed later on
            #    # Thomas: The solution is not to not integrate aloing the axis. Anyway there will not be any fluxes on the axis.
            #    weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(fluc.u),conj(self.X[0])),mean.u),self.n) )*self.x[1]*self.all_ds)
            #    weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
            #    weakForm.add(( -1j*fluc.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
                        
        else:
            ## ---- ALTERNATIVE: No integration by part, just one volume term
            logger.debug(" -> NOT using integration by parts for convection term.")
            # -- > Tensor implementation derived by hand
            weakForm += (-1j*iDot(iDot(iGrad(fluc.u), mean.rho*mean.u), iConj(X))).ufl_tens*J_hat*dx
            weakForm += (-1j*iDot(iDot(iGrad(mean.u), mean.rho*fluc.u), iConj(X))).ufl_tens*J_hat*dx
            weakForm += (-1j*iDot(iDot(iGrad(mean.u), fluc.rho*mean.u), iConj(X))).ufl_tens*J_hat*dx


        # ------------------------ Pressure gradient terms
        int_by_parts = True
        if int_by_parts:
            # Integrate pressure gradient boundary terms (resulting from integration by parts)
            weakForm += (1j*fluc.p*iDiv(iConj(X))).ufl_tens*J_hat*dx
            weakForm += (-1j*iDot(fluc.p*iConj(X), self.n)).ufl_tens*J_hat*self.all_ds
        else:
            # No integration by parts of the pressure term
            weakForm += -(1j*iDot(iGrad(fluc.p), iConj(X)) ).ufl_tens*J_hat*dx

        # ------------------------ Diffusion term
        # NOTE: In the current implementation of FELiCS, a mean.rhoean factor is missing
        #       in front of the viscosity. This error is kept for now for concistency,
        #       but it will need to be corrected. Thomas: The name of the variable is wrong, 
        #       the equations are correct. The nu is actually a mu. This needs to be corrected
        # NOTE: The diffusion term in the previous implementation of FELiCS neglects
        #       spatial gradients of the viscosity. => Sophie: This is corrected now, see
        #       'Fields/fieldProperties.py' for the definition of 'tau'.
        # NOTE: The boundary term from the integration by part is ignored. This should impose a 
        #       BC equivalent to stress-free BC
        

        weakForm += ( -1j*iInner(fluc.tau,iGrad(iConj(X)) )).ufl_tens*J_hat*dx
        #weakForm.add(( 1j*iDot(iDot(fluc.tau,self.n ),iConj(X))).ufl_tens*J_hat*self.all_ds)

        ## ---- Visc. 3: viscous BC terms for input-output analysis
        if self.param.Case.AnalysisMode in ['Input-Output']:
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:
                weakForm += ( 1j*iDot(iDot(fluc.tau,self.n),iConj(X)) ).ufl_tens*J_hat*self.ds(boundary_index)
                # Version with full viscous tensor (not assuming constant viscosity) --> Not working as expected for now
                #self.A_vf.add((1j*mean.nuTot*iDot(iDot(iGrad(fluc.u, self.m)+iT(iGrad(fluc.u, self.m)),),iConj(X))).ufl_tens*J_hat*self.ds(boundary_index))

    def addNonlinearExpression(self, weakForm, mean):
        """
        Add the nonlinear expression to the weak form.

        Parameters
        ----------
        weakForm : Form
            The weak form object to be updated.
        mean : MeanFields
            The mean fields object containing time-averaged variables.

        Notes
        -----
        - This method assembles the full nonlinear form of the momentum conservation equation. 
        - It supports integration by parts for convective and pressure terms and includes the diffusion term using the mean stress tensor.
        """

        X     = self.X
        J_hat = self.J_hat

        # ------------------------ Convective terms
        int_by_parts = True
        if int_by_parts and self.param.Case.CoordinateSystem=='Cartesian':
            logger.debug(" -> Using integration by parts for convection term.")
            # Volume term from integration by parts
            weakForm.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*mean.u)),mean.u) ).ufl_tens*J_hat*dx)
            # Boundary term from integration by parts
            weakForm.add(( -1j*mean.rho*iDot(iDot(iOuter(mean.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)
            #elif self.param.Case.CoordinateSystem =='Cylindrical':
            #    # In cyl , a singular term error arise for the boundary term in the tensor framework
            #    # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
            #    # This should be fixed later on
            #    # Thomas: The solution is not to not integrate aloing the axis. Anyway there will not be any fluxes on the axis.
            #    weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),mean.u),self.n) )*self.x[1]*self.all_ds)
                        
        else:
            ## ---- ALTERNATIVE: No integration by part, just one volume term
            logger.debug(" -> NOT using integration by parts for convection term.")
            # -- > Tensor implementation derived by hand
            weakForm.add(( -1j*iDot(iDot(iGrad(mean.u),mean.rho*mean.u),iConj(X)) ).ufl_tens*J_hat*dx)

        # ------------------------ Pressure gradient terms
        int_by_parts = True  
        if int_by_parts:
            # Integrate pressure gradient boundary terms (resulting from integration by parts)
            weakForm.add( (1j*mean.p*iDiv(iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add( (-1j*iDot(mean.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)

        else:
            # No integration by parts of the pressure term
            logger.error(' -> Pressure term without IbP not implemented in tensor framework.')
            raise Exception(' -> Pressure term without IbP not implemented in tensor framework.')

        # ------------------------ Diffusion term
        weakForm.add(( -1j*iInner(mean.tau,iGrad(iConj(X)) )).ufl_tens*J_hat*dx)


    def addBilinearExpression(self, weakForm, mean):
        """
        Add the bilinear convection term for incompressible flows.

        Parameters
        ----------
        weakForm : Form
            The weak form object to be updated.
        mean : MeanFields
            The mean fields object containing time-averaged variables.

        Notes
        -----
        - This method implements a simplified bilinear form for incompressible flow cases in the BOA project. 
        - Only the convective term is considered, using a strong formulation.
        """

        # Sophie: first and quick implementation for the BOA project. Only for incompressible flow and strong formulation (only convection term)
        u_bil = Tensor( mean._fieldDict['u_bilinear'], self.coordinateSystem)
        fluc  = self.fluc
        X     = self.X
        J_hat = self.J_hat

        weakForm.add(( -1j*iDot(iDot(iGrad(fluc.u),mean.rho*u_bil),iConj(X)) ).ufl_tens*J_hat*dx)
        weakForm.add(( -1j*iDot(iDot(iGrad(u_bil),mean.rho*fluc.u),iConj(X)) ).ufl_tens*J_hat*dx)


