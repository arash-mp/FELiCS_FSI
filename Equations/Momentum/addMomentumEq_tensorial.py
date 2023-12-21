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
    transpose
)
from tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT,
    iIdentity
)


from functions import printWarning, printError, printDebug

def addMomentumEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    momentum conservation equation, in tensorial framework.
    The current version is based on the addMomentumEq.py function
    that was previously defined using the index notation of ufl.
    For now all terms from the momentum equation have been grouped 
    into a single file instead of being in serparate files.

    NOTE:   Support for discontinuous Galerkin has been dropped
                
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(x)
    '''
    
    # Disclaimer
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')


    # ------------------------ Define the tensorial operators
    # Coordinate system
    coord = self.coord_sys
    
    ## Vector quantities
    #if fluc.u.ufl_shape[0] == 2:
    #    printDebug(param.debug, '--> Mom eq: fluctuations and mean flow 2D')
    #    x_tens = Tensor(as_vector((X[0], X[1], 0.0)), coord,containsTestFunction=True)
    #    u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord)
    #    u_m = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord)
    #
    #elif fluc.u.ufl_shape[0] == 3:
    #    printDebug(param.debug, '--> Mom eq: fluctuations and mean flow 3D')
    #    x_tens = Tensor(as_vector((X[0], X[1], X[2])), coord,containsTestFunction=True)
    #    u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord)
    #    u_m = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord)
    #    
    #else:
    #    printError('--> Mom eq: u has neither 2 or 3 dimensions: not implemented.') 
        
    x_tens = Tensor(X, coord,containsTestFunction=True)
    u_f = Tensor(fluc.u , coord)
    u_m = Tensor(mean.u, coord)
    # Scalar quantities
    rho_m = Tensor(mean.rho, coord)
    rho_f = Tensor(fluc.rho, coord)
    nutot_m = Tensor(mean.nuTot, coord)
    nulam_f = Tensor(fluc.nulam, coord)
    p_f = Tensor(fluc.p, coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)
        

    # ------------------------ Time derivative term
    self.B_vf.add(( iDot(rho_m*u_f,iConj(x_tens)) ).ufl_tens*coord.J_hat*dx)


    # ------------------------ Convective terms
    int_by_parts = True
    if int_by_parts:
        # Volume term from integration by parts
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(x_tens),rho_m*u_m),self.m),u_f) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(x_tens),rho_m*u_f),0),u_m) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*iDot(iDiv(iOuter(iConj(x_tens),rho_f*u_m),0),u_m) ).ufl_tens*coord.J_hat*dx)
        # Boundary term from integration by parts
        if param.Case.CoordinateSystem =='Cartesian':
            print( rho_m*iDot(iDot(iOuter(u_f,iConj(x_tens)),u_m),nbc_tens) )
            self.A_vf.add(( -1j*rho_m*iDot(iDot(iOuter(u_f,iConj(x_tens)),u_m),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
            self.A_vf.add(( -1j*rho_m*iDot(iDot(iOuter(u_m,iConj(x_tens)),u_f),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
            self.A_vf.add(( -1j*rho_f*iDot(iDot(iOuter(u_m,iConj(x_tens)),u_m),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
        elif param.Case.CoordinateSystem =='Cylindrical':
            # In cyl coords, a singular term error arise for the boundary term in the tensor framework
            # Because there is no Nabla operator in the boundary term we can use the ufl operator and avoid this error
            # This should be fixed later on
            self.A_vf.add(( -1j*mean.rho*dot(dot(outer(conj(fluc.u),conj(X)),mean.u),as_vector((self.n_BC[0],self.n_BC[1],0.0))) )*self.x[1]*self.all_ds)
            self.A_vf.add(( -1j*mean.rho*dot(dot(outer(conj(mean.u),conj(X)),fluc.u),as_vector((self.n_BC[0],self.n_BC[1],0.0))) )*self.x[1]*self.all_ds)
            self.A_vf.add(( -1j*fluc.rho*dot(dot(outer(conj(mean.u),conj(X)),fluc.u),as_vector((self.n_BC[0],self.n_BC[1],0.0))) )*self.x[1]*self.all_ds)
        else:
            printError('Coord. syst not yet implemented in tensor framework.')
                    
    else:
        ## ---- ALTERNATIVE: No integration by part, just one volume term
        printDebug(param.debug, '--> Mom eq: convection term NOT integrated by part')
        # -- > Tensor implementation derived by hand
        self.A_vf.add(( -1j*iDot(iDot(iGrad(u_f,self.m),rho_m*u_m),iConj(x_tens)) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( -1j*iDot(iDot(iGrad(u_m),rho_m*u_f),iConj(x_tens)) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( -1j*iDot(iDot(iGrad(u_m),rho_f*u_m),iConj(x_tens)) ).ufl_tens*coord.J_hat*dx)


    # ------------------------ Pressure gradient terms
    int_by_parts = True  
    if int_by_parts:
        # Integrate pressure gradient boundary terms (resulting from integration by parts)
        self.A_vf.add( (1j*p_f*iDiv(iConj(x_tens),self.m) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add( (-1j*iDot(p_f*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)

    else:
        # No integration by parts of the pressure term
        printError('--> Mom eq: Pressure term without IbP not implemented in tensor framework.')
        # -- > Previous implementation
        #self.A_vf.add(1j * -self.R*inner(grad(fluc.p), X)*dx)
        # Additionnal terms for m > 0
        #if not param.Case.m == 0: # Should this not be ony for cyl. coords.??
        #    self.A_vf.add(+conj(X)[2]*self.m*fluc.p*dx)



    # ------------------------ Diffusion terms
    # NOTE: In the current implementation of FELiCS, a rho_mean factor is missing
    # in front of the viscosity. This error is kept for now for concistency,
    # but it will need to be corrected.
    # NOTE: The diffusion term in the previous implementation of FELiCS neglects
    # spatial gradients of the viscosity.
    # NOTE: The boundary term from the integration by part is ignored. This should impose a 
    # BC equivalent to stress-free BC
    
    
    ## ---- Visc. 1: viscous terms for incompressible flow
    # Term corresponding to:     div(mu*gradT(u)) = div(mu_mean*grad(u_fluc)) + div(mu_fluc*grad(u_mean))
    self.A_vf.add(( -1j*nutot_m*iInner(iGrad(u_f,self.m),iGrad(iConj(x_tens),self.m) )).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( -1j*nulam_f*iInner(iGrad(u_m),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
    
    # In case we decide not to neglect the gradient of nutot (still need debugging):
    # Term corresponding to:    div(mu*grad^T(u)) = div(mu_mean*grad^T(u_fluc)) + div(mu_fluc*grad^T(u_mean))
    # self.A_vf.add(( -1j*nutot_m*iInner(iT(iGrad(u_f,self.m)),iGrad(iConj(x_tens),-self.m)) ).ufl_tens*coord.J_hat*dx)
    # self.A_vf.add(( -1j*nulam_f*iInner(iT(iGrad(u_m)),iGrad(iConj(x_tens),-self.m)) ).ufl_tens*coord.J_hat*dx)
    # Do we need to add the following boundary term then? --> Singular value error!
    # self.A_vf.add(( 1j*nutot_m*iDot(iDot(iT(iGrad(u_f,self.m)),nbc_tens),iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
    # self.A_vf.add(( 1j*nulam_f*iDot(iDot(iT(iGrad(u_m)),nbc_tens),iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)


    ## ---- Visc. 2: viscous terms for compressible flow
    if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
        # The viscous diffusion is still integrated by parts and the resulting boundary 
        # terms are neglected to impose a Neumann BC
        # Adding term corresponding to: div(mu*grad^T(u)) = div(mu_mean*grad^T(u_fluc)) + div(mu_fluc*grad^T(u_mean))
        self.A_vf.add(( -1j*nutot_m*iInner(iT(iGrad(u_f,self.m)),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( -1j*nulam_f*iInner(iT(iGrad(u_m)),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
        # Adding terms corresponding to: -2/3*div(mu*(div(u)*I) = -2/3*div(mu_mean*(div(u_fluc)*I) -2/3*div(mu_fluc*(div(u_mean)*I)
        self.A_vf.add(( 1j*2.0/3.0*nutot_m*iInner(iDiv(u_f,self.m)*iIdentity(iGrad(u_f)),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*2.0/3.0*nulam_f*iInner(iDiv(u_m)*iIdentity(iGrad(u_f)),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
    

    ## ---- Visc. 3: viscous BC terms for input-output analysis
    if param.Case.AnalysisMode in ['Input-Output']:
        for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
            self.A_vf.add(( 1j*nutot_m*iDot(iDot(iGrad(u_f,self.m),nbc_tens),iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(boundary_index))
            # Version with full viscous tensor (not assuming constant viscosity) --> Not working as expected for now
            #self.A_vf.add((1j*nutot_m*iDot(iDot(iGrad(u_f, self.m)+iT(iGrad(u_f, self.m)),nbc_tens),iConj(x_tens))).ufl_tens*coord.J_hat*self.ds(boundary_index))
            
