from ufl import (
    dx,
    conj,
    Identity,
    i,
    j,
    k,
    Dx,
    as_tensor,
    inner,
    grad,
    dot,
    outer,
    transpose,
    Constant,
)
from FELiCS.Misc.tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT,
    iIdentity,
)


from FELiCS.Misc.functions import printWarning, printError, printDebug

from .EquationTemplate import EquationTemplate

class MomentumEquation(EquationTemplate):

    def __init__(self,eqColl,fluc,X,param):
 
        # Disclaimer
        if param.NumericalScheme in ['Discontinuous Galerkin']:
            printError('Discontinuous Galerkin not implemented in tensorial framework.')
    
        # initialize variables in template class
        super().__init__(eqColl,fluc,X,param)


    def addWeightMatrixExpression(self,weakForm,mean):
        # Time derivative term
        weakForm.add(( iDot(mean.rho*self.fluc.u,iConj(self.X)) ).ufl_tens*self.J_hat*dx)

    def addLinearExpression(self,weakForm,mean):
        '''
        This function builds the weak form of the linearized
        momentum conservation equation, in tensorial framework.
        The current version is based on the addMomentumEq.py function
        that was previously defined using the index notation of ufl.
        For now all terms from the momentum equation have been grouped 
        into a single file instead of being in serparate files.
        '''
      
        fluc  = self.fluc
        X     = self.X
        J_hat = self.J_hat

        # ------------------------ Convective terms
        int_by_parts = False
        if int_by_parts:
            # Volume term from integration by parts
            weadForm.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*mean.u)),fluc.u) ).ufl_tens*J_hat*dx)
            weadForm.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*fluc.u)),mean.u) ).ufl_tens*J_hat*dx)
            weadForm.add(( 1j*iDot(iDiv(iOuter(iConj(X),fluc.rho*mean.u)),mean.u) ).ufl_tens*J_hat*dx)
            # Boundary term from integration by parts
            if self.param.Case.CoordinateSystem =='Cartesian':
                weakForm.add(( -1j*mean.rho*iDot(iDot(iOuter(fluc.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)
                weakForm.add(( -1j*mean.rho*iDot(iDot(iOuter(mean.u,iConj(X)),fluc.u),self.n) ).ufl_tens*J_hat*self.all_ds)
                weakForm.add(( -1j*fluc.rho*iDot(iDot(iOuter(mean.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)
            elif self.param.Case.CoordinateSystem =='Cylindrical':
                # In cyl , a singular term error arise for the boundary term in the tensor framework
                # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
                # This should be fixed later on
                # Thomas: The solution is not to not integrate aloing the axis. Anyway there will not be any fluxes on the axis.
                weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(fluc.u),conj(self.X[0])),mean.u),self.n) )*self.x[1]*self.all_ds)
                weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
                weakForm.add(( -1j*fluc.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
            else:
                printError('Coord. syst not yet implemented in tensor framework.')
                        
        else:
            ## ---- ALTERNATIVE: No integration by part, just one volume term
            printDebug(self.param.debug, '-- -> Mom eq: convection term NOT integrated by part')
            # -- > Tensor implementation derived by hand
            weakForm.add(( -1j*iDot(iDot(iGrad(fluc.u),mean.rho*mean.u),iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( -1j*iDot(iDot(iGrad(mean.u),mean.rho*fluc.u),iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( -1j*iDot(iDot(iGrad(mean.u),fluc.rho*mean.u),iConj(X)) ).ufl_tens*J_hat*dx)


        # ------------------------ Pressure gradient terms
        int_by_parts = True  
        if int_by_parts:
            # Integrate pressure gradient boundary terms (resulting from integration by parts)
            weakForm.add( (1j*fluc.p*iDiv(iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add( (-1j*iDot(fluc.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)

        else:
            # No integration by parts of the pressure term
            printError('-- -> Mom eq: Pressure term without IbP not implemented in tensor framework.')



        # ------------------------ Diffusion term
        # NOTE: In the current implementation of FELiCS, a mean.rhoean factor is missing
        # in front of the viscosity. This error is kept for now for concistency,
        # but it will need to be corrected. Thomas: The name of the variable is wrong, the equations are correct. The nu is actually a mu. This needs to be corrected
        # NOTE: The diffusion term in the previous implementation of FELiCS neglects
        # spatial gradients of the viscosity. Thomas: It does not!!! Should be correct
        # NOTE: The boundary term from the integration by part is ignored. This should impose a 
        # BC equivalent to stress-free BC
        
        weakForm.add(( -1j*iInner(fluc.tau,iGrad(iConj(X)) )).ufl_tens*J_hat*dx)

        ## ---- Visc. 3: viscous BC terms for input-output analysis
        if self.param.Case.AnalysisMode in ['Input-Output']:
            for boundary_index in self.param.IOResolvent.ForcingBoundaryIndices:
                weakForm.add(( 1j*mean.nuTot*iDot(iDot(iGrad(fluc.u),self.n),iConj(X)) ).ufl_tens*J_hat*self.ds(boundary_index))
                # Version with full viscous tensor (not assuming constant viscosity) --> Not working as expected for now
                #self.A_vf.add((1j*mean.nuTot*iDot(iDot(iGrad(fluc.u, self.m)+iT(iGrad(fluc.u, self.m)),),iConj(X))).ufl_tens*J_hat*self.ds(boundary_index))


    def addNonlinearExpression(self, weakForm, mean):
      
        X     = self.X
        J_hat = self.J_hat

        # ------------------------ Convective terms
        int_by_parts = False
        if int_by_parts:
            # Volume term from integration by parts
            weakForm.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*mean.u)),mean.u) ).ufl_tens*J_hat*dx)
            # Boundary term from integration by parts
            if self.param.Case.CoordinateSystem =='Cartesian':
                weakForm.add(( -1j*mean.rho*iDot(iDot(iOuter(mean.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)
            elif self.param.Case.CoordinateSystem =='Cylindrical':
                # In cyl , a singular term error arise for the boundary term in the tensor framework
                # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
                # This should be fixed later on
                # Thomas: The solution is not to not integrate aloing the axis. Anyway there will not be any fluxes on the axis.
                weakForm.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),mean.u),self.n) )*self.x[1]*self.all_ds)
            else:
                printError('Coord. syst not yet implemented in tensor framework.')
                        
        else:
            ## ---- ALTERNATIVE: No integration by part, just one volume term
            printDebug(self.param.debug, '-- -> Mom eq: convection term NOT integrated by part')
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
            printError('-- -> Mom eq: Pressure term without IbP not implemented in tensor framework.')

        # ------------------------ Diffusion term
        # NOTE: In the current implementation of FELiCS, a mean.rhoean factor is missing
        # in front of the viscosity. This error is kept for now for concistency,
        # but it will need to be corrected. Thomas: The name of the variable is wrong, the equations are correct. The nu is actually a mu. This needs to be corrected
        # NOTE: The diffusion term in the previous implementation of FELiCS neglects
        # spatial gradients of the viscosity. Thomas: It does not!!! Should be correct
        # NOTE: The boundary term from the integration by part is ignored. This should impose a 
        # BC equivalent to stress-free BC
        weakForm.add(( -1j*iInner(mean.tau,iGrad(iConj(X)) )).ufl_tens*J_hat*dx)


    def addBilinearExpression(self, weakForm, mean):
        # Sophie: first and quick implementation for the BOA project. Only for incompressible flow and strong formulation (only convection term)
        u_bil = Tensor( mean._fieldDict['u_bilinear'], self.coordinateSystem)
        fluc  = self.fluc
        X     = self.X
        J_hat = self.J_hat

        weakForm.add(( -1j*iDot(iDot(iGrad(fluc.u),mean.rho*u_bil),iConj(X)) ).ufl_tens*J_hat*dx)
        weakForm.add(( -1j*iDot(iDot(iGrad(u_bil),mean.rho*fluc.u),iConj(X)) ).ufl_tens*J_hat*dx)


