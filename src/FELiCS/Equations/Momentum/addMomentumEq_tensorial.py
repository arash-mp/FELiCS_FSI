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
from FELiCS.tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iDotT,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT,
    iIdentity,
)


from FELiCS.functions import printWarning, printError, printDebug
from FELiCS.Equations.Mass.addMassEq_tensorial import addMassEq

def addMomentumEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    momentum conservation equation, in tensorial framework.
    The current version is based on the addMomentumEq.py function
    that was previously defined using the index notation of ufl.
    For now all terms from the momentum equation have been grouped 
    into a single file instead of being in serparate files.
    '''
    
    # Disclaimer
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')


    # ------------------------ Define the tensorial operators
    J_hat = self._coordinateSystem.J_hat

    # ------------------------ Time derivative term
    self.B_vf.add(( iDot(mean.rho*fluc.u,iConj(X)) ).ufl_tens*J_hat*dx)
    # ------------------------ Convective terms
    int_by_parts = True
    if param.Case.CoordinateSystem =='Cartesian' and int_by_parts:
        # Volume term from integration by parts

        # version0: standard formulation 
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*mean.u)),fluc.u) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(X),mean.rho*fluc.u)),mean.u) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(X),fluc.rho*mean.u)),mean.u) ).ufl_tens*J_hat*dx)

        # version1: Done by Sophie: using conservative form for stabilization
        #self.A_vf.add(( 1j*iDot(iDotT(iGrad(iConj(X)),mean.rho*mean.u ), fluc.u) ).ufl_tens*J_hat*dx)
        #self.A_vf.add(( 1j*iDot(iDotT(iGrad(iConj(X)),mean.rho*fluc.u ), mean.u) ).ufl_tens*J_hat*dx)
        #self.A_vf.add(( 1j*iDot(iDotT(iGrad(iConj(X)),fluc.rho*mean.u ), mean.u) ).ufl_tens*J_hat*dx)

        # Boundary term from integration by parts
        self.A_vf.add(( -1j*mean.rho*iDot(iDot(iOuter(fluc.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)
        self.A_vf.add(( -1j*mean.rho*iDot(iDot(iOuter(mean.u,iConj(X)),fluc.u),self.n) ).ufl_tens*J_hat*self.all_ds)
        self.A_vf.add(( -1j*fluc.rho*iDot(iDot(iOuter(mean.u,iConj(X)),mean.u),self.n) ).ufl_tens*J_hat*self.all_ds)

        #elif param.Case.CoordinateSystem =='Cylindrical':
        #    # In cyl , a singular term error arise for the boundary term in the tensor framework
        #    # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
        #    # This should be fixed later on
        #    # Thomas: The solution is not to not integrate aloing the axis. Anyway there will not be any fluxes on the axis.
        #    self.A_vf.add(( -1j*mean.rho*dot(dot(outer(conj(fluc.u),conj(self.X[0])),mean.u),self.n) )*self.x[1]*self.all_ds)
        #    self.A_vf.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
        #    self.A_vf.add(( -1j*fluc.rho*dot(dot(outer(conj(mean.u),conj(self.X[0])),fluc.u),self.n) )*self.x[1]*self.all_ds)
        #else:
        #    printError('Coord. syst not yet implemented in tensor framework.')
                    
    elif param.Case.CoordinateSystem == 'Cylindrical' or not int_by_parts:
        ## ---- ALTERNATIVE: No integration by part, just one volume term
        printDebug(param.debug, '-- -> Mom eq: convection term NOT integrated by part')
        # -- > Tensor implementation derived by hand
        # Done by Sophie: using conservative form for stabilization
        #self.A_vf.add(( -1j * iDot( iDiv( iOuter(mean.rho*mean.u, fluc.u)), iConj(X)) ).ufl_tens*J_hat*dx)
        #self.A_vf.add(( -1j * iDot( iDiv( iOuter(mean.rho*fluc.u, mean.u)), iConj(X)) ).ufl_tens*J_hat*dx)
        #self.A_vf.add(( -1j * iDot( iDiv( iOuter(fluc.rho*mean.u, mean.u)), iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( -1j*iDot(iDot(iGrad(fluc.u),mean.rho*mean.u),iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( -1j*iDot(iDot(iGrad(mean.u),mean.rho*fluc.u),iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( -1j*iDot(iDot(iGrad(mean.u),fluc.rho*mean.u),iConj(X)) ).ufl_tens*J_hat*dx)

    else:
        printError('Coord. syst not yet implemented in tensor framework.')

    # ------------------------ Pressure gradient terms
    int_by_parts = True  
    if int_by_parts:
        # Integrate pressure gradient boundary terms (resulting from integration by parts)
        self.A_vf.add( (1j*fluc.p*iDiv(iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add( (-1j*iDot(fluc.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)

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
    
    self.A_vf.add(( -1j*iInner(fluc.tau,iGrad(iConj(X)) )).ufl_tens*J_hat*dx)

    ## ---- Visc. 3: viscous BC terms for input-output analysis
    if param.Case.AnalysisMode in ['Input-Output']:
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            self.A_vf.add(( 1j*mean.nuTot*iDot(iDot(iGrad(fluc.u),self.n),iConj(X)) ).ufl_tens*J_hat*self.ds(boundary_index))
            # Version with full viscous tensor (not assuming constant viscosity) --> Not working as expected for now
            #self.A_vf.add((1j*mean.nuTot*iDot(iDot(iGrad(fluc.u, self.m)+iT(iGrad(fluc.u, self.m)),),iConj(X))).ufl_tens*J_hat*self.ds(boundary_index))
           



