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
    outer
)
from tensor_utils import (
    Tensor,
    CoordinateSystem,
    as_vector,
    as_matrix,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT
)
from functions import printWarning, printError

def addMomentumEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    momentum conservation equation, in tensorial framework.
    The current version is based on the addMomentumEq.py function
    that was previously defined using the index notation of ufl.
    For now all terms from the momentum equation have been grouped 
    into a single file instead of being in serparate files.

    NOTE:
        (i)     This version doesn't work in cylindrical coordinates
        (ii)    It is an early attempt at integrating the tensor
                framework in FELiCS, so bugs/errors are expected.
        (iii)   Support for discontinuous Galerkin has been dropped.
        (iv)    Viscosity fluctuations are currently neglected.
        (v)     This version is only for an INCOMPRESSIBLE flow, 
                and does not include density fluctuations.
        (vi)    The "FLAG_TENS" is used to activate the tensor
                framework in specific terms of the eq. It is possible
                to combine tensorial terms with index-notation ones
    '''

    # ------------------------ Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED "mesh_dims" and "as_vector(x1, x2, 0.0)" to cartesian coords
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        x_tens = Tensor(as_vector((X[0], X[1], 0.0)), coord_sys)
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # HARDCODED FOR CYL COORD IN FELiCS, DEFINED AS (Z, R, PHI)!
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        x_tens = Tensor(as_vector((X[0], X[1], X[2])), coord_sys)
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord_sys)
        
    else:
        printError('Coord. syst not yet implemented in tensor framework.')
        
    rho_mean_tens = Tensor(mean.rho, coord_sys)
    nutot_mean_tens = Tensor(mean.nuTot, coord_sys)
    nulam_fluc_tens = Tensor(fluc.nulam, coord_sys)
    rho_fluc_tens = Tensor(fluc.rho, coord_sys)
    p_fluc_tens = Tensor(fluc.p, coord_sys)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')
        

    # ------------------------ Time derivative term
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (not working for cyl coord)
        self.B_vf.add( (iDot(rho_mean_tens*u_fluc_tens, iConj(x_tens))).ufl_tens * coord_sys.J_hat * dx)
        
    else:
        # -- > Previous implementation
        self.B_vf.add(mean.rho*self.R*fluc.u[i]*conj(X)[i] *dx)


    # ------------------------ Convective terms
    FLAG_TENS = True
    int_by_parts = True
    
    if int_by_parts:
        
        if FLAG_TENS:
            # OPT1. -- > Tensor implementation (derived by hand by Jens)
            # -- Volume term from integration by parts
            self.A_vf.add((1j*iDot(iDiv(iOuter(iConj(x_tens),rho_mean_tens*u_mean_tens), -self.m),u_fluc_tens)).ufl_tens*coord_sys.J_hat*dx)
            self.A_vf.add((1j*iDot(iDiv(iOuter(iConj(x_tens),rho_mean_tens*u_fluc_tens), 0),u_mean_tens)).ufl_tens*coord_sys.J_hat*dx)
            self.A_vf.add((1j*iDot(iDiv(iOuter(iConj(x_tens),rho_fluc_tens*u_mean_tens), 0),u_mean_tens)).ufl_tens*coord_sys.J_hat*dx)
            
            # -- Boundary term from integration by parts
            if param.Case.CoordinateSystem =='Cartesian':
                # We might need to use the ufl form in cartesian coord with m > 0 (not tested yet)
                self.A_vf.add( (-1j * rho_mean_tens*iDot(iDot( iOuter(u_fluc_tens, iConj(x_tens)), u_mean_tens), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
                self.A_vf.add( (-1j * rho_mean_tens*iDot(iDot( iOuter(u_mean_tens, iConj(x_tens)), u_fluc_tens), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
                self.A_vf.add( (-1j * rho_fluc_tens*iDot(iDot( iOuter(u_mean_tens, iConj(x_tens)), u_mean_tens), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
                
            elif param.Case.CoordinateSystem =='Cylindrical':
                # In cyl coords, a singular term error arise for the above term in the tensor framework
                # Because there is no Nabla in this term we can use the ufl operators, which seems this error
                self.A_vf.add( (-1j * mean.rho*dot(dot( outer(conj(fluc.u), conj(X)), mean.u), as_vector((self.n_BC[0], self.n_BC[1], 0.0)) )) *self.x[1]*self.all_ds)
                self.A_vf.add( (-1j * mean.rho*dot(dot( outer(conj(mean.u), conj(X)), fluc.u), as_vector((self.n_BC[0], self.n_BC[1], 0.0)) )) *self.x[1]*self.all_ds)
                self.A_vf.add( (-1j * fluc.rho*dot(dot( outer(conj(mean.u), conj(X)), fluc.u), as_vector((self.n_BC[0], self.n_BC[1], 0.0)) )) *self.x[1]*self.all_ds)
                
            else:
                printError('Coord. syst not yet implemented in tensor framework.')
        
        else:
            # -- > Previous implementation
            # -- Volume term from integration by parts
            I = Identity(fluc.u.geometric_dimension())
            self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * mean.u[j] ,k) * fluc.u[i] * I[k,j] * dx)
            self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * fluc.u[j] ,k) * mean.u[i] * I[k,j] * dx)
            self.A_vf.add(1j * Dx(self.R * conj(X[i]) * fluc.rho * mean.u[j] ,k) * mean.u[i] * I[k,j] * dx)
            
            # -- Boundary term from integration by parts
            self.A_vf.add(1j * (- self.R * conj(X[i]) * mean.rho * mean.u[j] * fluc.u[i] * self.n_BC[k] * I[k,j]\
			    - self.R * conj(X[i]) * mean.rho * fluc.u[j] * mean.u[i] * self.n_BC[k] * I[k,j]\
			    - self.R * conj(X[i]) * fluc.rho * mean.u[j] * mean.u[i] * self.n_BC[k] * I[k,j])*self.all_ds)
            
            # -- Additionnal terms for cyl. coord
            if param.Case.CoordinateSystem in ['Cylindrical']:
                self.A_vf.add(1j * - conj(X)[2] * mean.rho *mean.u[1] * fluc.u[2] * dx)
                self.A_vf.add(1j * - conj(X)[2] * mean.rho *fluc.u[1] * mean.u[2] * dx)
                self.A_vf.add(1j * - conj(X)[2] * fluc.rho *mean.u[1] * mean.u[2] * dx)

                self.A_vf.add(1j * conj(X)[1] * mean.u[2] * fluc.u[2] * mean.rho * dx\
					  + 1j * conj(X)[1] * fluc.u[2] * mean.u[2] * mean.rho * dx\
					  + 1j * conj(X)[1] * mean.u[2] * mean.u[2] * fluc.rho * dx)

                # -- Additionnal terms for cyl. coord and m > 0
                if not param.Case.m == 0:
                    self.A_vf.add(mean.rho * conj(X)[i] * param.Case.m * mean.u[2] * fluc.u[i] * dx)
                    
    else:
         ## ---- ALTERNATIVE: No integration by part, just one volume term
        if not param.Case.m == 0:
            printError('--> Mom eq: Convection term without IbP not implemented for m > 0.')     
             
        if FLAG_TENS:
            # -- > Tensor implementation derived by hand
            self.A_vf.add( (-1j*iDot( iDot(iGrad(u_fluc_tens, self.m), rho_mean_tens*u_mean_tens), iConj(x_tens)) ).ufl_tens*coord_sys.J_hat* dx)
            self.A_vf.add( (-1j*iDot( iDot(iGrad(u_mean_tens), rho_mean_tens*u_fluc_tens), iConj(x_tens)) ).ufl_tens*coord_sys.J_hat* dx)
            self.A_vf.add( (-1j*iDot( iDot(iGrad(u_mean_tens), rho_fluc_tens*u_mean_tens), iConj(x_tens)) ).ufl_tens*coord_sys.J_hat* dx)
            
        else:
            # Index notation implementation (m > 0 terms missing!)
            I = Identity( fluc.u.geometric_dimension() )
            self.A_vf.add(-1j * self.R*conj(X[i])*mean.rho*mean.u[j]*Dx(fluc.u[i], k)*I[k,j] * dx)
            self.A_vf.add(-1j * self.R*conj(X[i])*mean.rho*fluc.u[j]*Dx(mean.u[i], k)*I[k,j] * dx)
            self.A_vf.add(-1j * self.R*conj(X[i])*fluc.rho*mean.u[j]*Dx(mean.u[i], k)*I[k,j] * dx)
            if param.Case.CoordinateSystem in ['Cylindrical']:
                self.A_vf.add(1j * conj(X)[1]*mean.rho*mean.u[2]*fluc.u[2] * dx)
                self.A_vf.add(1j * -conj(X)[2]*mean.rho*mean.u[1]*fluc.u[2] * dx)
                
                self.A_vf.add(1j * conj(X)[1]*mean.rho*fluc.u[2]*mean.u[2] * dx)
                self.A_vf.add(1j * -conj(X)[2]*mean.rho*fluc.u[1]*mean.u[2] * dx)
                
                self.A_vf.add(1j * conj(X)[1]*fluc.rho*mean.u[2]*mean.u[2] * dx)
                self.A_vf.add(1j * -conj(X)[2]*fluc.rho*mean.u[1]*mean.u[2] * dx)


    # ------------------------ Pressure gradient terms
    FLAG_TENS = True
    int_by_parts = True
    
    if int_by_parts:
        # Integrate pressure gradient boundary terms (resulting from integration by parts)
        if FLAG_TENS:
            # -- > Tensor implementation
            self.A_vf.add( (1j*p_fluc_tens*iDiv(iConj(x_tens), -self.m)).ufl_tens*coord_sys.J_hat*dx)
            self.A_vf.add( (-1j*iDot(p_fluc_tens*iConj(x_tens), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
            
        else:
            # -- > Previous implementation
            I = Identity( fluc.u.geometric_dimension() )
            self.A_vf.add(1j * Dx(conj(X[i])*self.R,j)*fluc.p * I[i,j] * dx)
            self.A_vf.add(1j * -self.R*fluc.p*self.n_BC[j]*conj(X[i])* I[i,j] * self.all_ds)
            
            # -- Additionnal terms for m > 0
            if not param.Case.m == 0: # Should this not be ony for cyl. coords.??
                self.A_vf.add(+conj(X)[2]*self.m*fluc.p*dx)

    else:
        # No integration by parts of the pressure term
        if FLAG_TENS:
            printError('--> Mom eq: Pressure term without IbP not implemented in tensor framework.')
            
        else:
            # -- > Previous implementation
            self.A_vf.add(1j * -self.R*inner(grad(fluc.p), X)*dx)
            
            # -- Additionnal terms for m > 0
            if not param.Case.m == 0: # Should this not be ony for cyl. coords.??
                self.A_vf.add(+conj(X)[2]*self.m*fluc.p*dx)



    # ------------------------ Diffusion terms
    # NOTE: In the current implementation of FELiCS, a rho_mean factor is missing
    # in front of the viscosity. This error is kept for now for concistency,
    # but it will need to be corrected.
    # NOTE: The diffusion term in the previous implementation of FELiCS neglects
    # spatial gradients of the viscosity.
    # NOTE: The boundary term from the integration by part is ignored. This should impose a 
    # BC equivalent to stress-free BC
    
    ## ---- Visc. 1: viscous terms for incompressible flow
    FLAG_TENS = True
    
    if FLAG_TENS:
        # -- > Tensor implementation
        self.A_vf.add((-1j*nutot_mean_tens*iInner(iGrad(u_fluc_tens,self.m), iGrad(iConj(x_tens),-self.m))).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add((-1j*nulam_fluc_tens*iInner(iGrad(u_mean_tens), iGrad(iConj(x_tens),-self.m))).ufl_tens*coord_sys.J_hat*dx)
        # NOTE: The term with nulam_fluc has not been validated yet
        
    else:
        # -- > Previous implementation
        # I,III: div(mean.mu*grad(fluc.u)) in 2D after Ibp (boundary integral neglected): 
        self.A_vf.add(1j * -self.R*mean.nuTot*Dx(fluc.u[i],j)*Dx(conj(X)[i],j)*dx )
        # div(fluc.mu*grad(mean.u)) in 2D after Ibp (integral of boundary terms neglected):
        self.A_vf.add(1j * -fluc.nulam*mean.u[i].dx(j)*conj(X)[i].dx(j)*dx )
        
        # Additional terms for flow in cyl. coord
        if param.Case.CoordinateSystem =='Cylindrical':
            self.A_vf.add(1j * -mean.nuTot/self.R*conj(X)[2]*fluc.u[2]*dx)
            self.A_vf.add(1j * -mean.nuTot.dx(1)*conj(X)[2]*fluc.u[2]*dx)
            self.A_vf.add(1j * -mean.nuTot*conj(X)[1]/self.R*fluc.u[1]*dx)
            
            # Additional terms for m > 0
            if not param.Case.m == 0:
                self.A_vf.add(-2*mean.nuTot/self.R*self.m*conj(X)[2]*fluc.u[1]*dx)#
                self.A_vf.add(-conj(X)[2]*mean.nuTot.dx(0)*self.m*fluc.u[0]*dx)
                #self.A_vf.add(1j * -X*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)            #III  -> A_real ??? (CA)
                self.A_vf.add(-conj(X)[2]*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)
                self.A_vf.add(1j * -mean.nuTot*self.m**2/self.R*inner(fluc.u, X)*dx)#
                #self.A_vf.add(1j * -mean.nuTot*self.m**2/self.R*conj(X)[i]*fluc.u[i]*dx)# this line replaces the above which uses inner
                self.A_vf.add(2*mean.nuTot*self.m/self.R*conj(X)[1]*fluc.u[2]*dx)#

    
    ## ---- Visc. 2: viscous terms for compressible flow
    if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
        if FLAG_TENS:
            printWarning('--> Mom eq: Compressible mom. eq. not implemented in tensor framework. Currently relies on index notation.')

        # -- > Previous implementation
        I = Identity( fluc.u.geometric_dimension() )
        # II, III: div(mu*grad^T(fluc.u)) in 2D after Ibp (boundary integral neglected):
        self.A_vf.add(1j * -self.R*mean.nuTot*Dx(fluc.u[j],i)*Dx(conj(X)[i],j)*dx )
        
        # IV: div(mean.mu*(div(u)*I)
        self.A_vf.add(1j * mean.nuTot*2.0/3.0*Dx(fluc.u[k],k)*Dx(self.R*conj(X)[i],j)*I[i,j]*dx ) # Change CA:
        #div(fluc.mu*grad^T(mean.u)) in 2D after Ibp (integral of boundary terms neglected):#
        #CAUTION: NOT IMPLEMENTED IN CYLINDRICAL COORDINATES!
        self.A_vf.add(1j * -fluc.nulam*mean.u[j].dx(i)*conj(X)[i].dx(j)*dx )
        self.A_vf.add(1j * 2.0/3.0 * fluc.nulam*mean.u[k].dx(k)*I[i,j]*conj(X)[i].dx(j)*dx )
        
        # Additional terms in cylindrical coordinates
        if param.Case.CoordinateSystem =='Cylindrical':
            self.A_vf.add(1j * -mean.nuTot/self.R*fluc.u[1]*conj(X)[1]*dx)              #II(1)
            self.A_vf.add(1j * -2.0/3.0*fluc.u[1]*Dx(mean.nuTot,i)*conj(X)[i]*dx)       #IV(1)
            self.A_vf.add(1j * -2.0/3.0*mean.nuTot*Dx(fluc.u[1],i)*conj(X)[i]*dx)       #IV(2)
            self.A_vf.add(1j * 2.0/3.0*mean.nuTot/self.R*fluc.u[1]*conj(X)[1]*dx)
            
            # Extra terms for m > 0
            if not param.Case.m == 0:
                self.A_vf.add(-self.m*mean.nuTot*Dx(fluc.u[2],i)*conj(X)[i]*dx)		    #II(2)
                self.A_vf.add(self.m/self.R*mean.nuTot*fluc.u[2]*conj(X)[1]*dx)		    #II(3)
                
                self.A_vf.add(2.0/3.0*self.m*fluc.u[2]*Dx(mean.nuTot,i)*conj(X)[i]*dx)  #IV(4)
                self.A_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[2],i)*conj(X)[i]*dx)  #IV(5)
                self.A_vf.add(-2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[2]*conj(X)[1]*dx)	
    

    ## ---- Visc. 3: viscous BC terms for input-output analysis
    FLAG_TENS = True
    if param.Case.AnalysisMode in ['Input-Output']:
        if FLAG_TENS:
            for boundary_index in param.IOResolvent.ForcingBoundaryIndices:
                self.A_vf.add((1j*nutot_mean_tens*iDot(iDot(iGrad(u_fluc_tens, self.m),nbc_tens),iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.ds(boundary_index))
            
        else:
            # -- > Previous implementation
            for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

                # Add physical boundary terms in imaginary part, where the forcing is applied...
                self.A_vf.add(1j * -self.R*mean.nuTot*
                            (-conj(X)[0] * (self.n_BC[0] * (fluc.u[0].dx(0)) +
                                            self.n_BC[1] * (fluc.u[0].dx(1))) -
                                conj(X)[1] *
                                            (self.n_BC[0] * (fluc.u[1].dx(0)) +
                                            self.n_BC[1] * (fluc.u[1].dx(1))) )
                                * self.ds(boundary_index))

                # Add stabilization terms on imaginary part according
                # to Baumann and Oden JFM 2016 vol 798
                # self.A_vf.add(1j * -self.R*mean.nuTot*((self.n_BC[0] * (conj(X)[0].dx(0)) +
                #                         self.n_BC[1] *   (conj(X)[0].dx(1))) * (fluc.u[0]-mean.u_forcing_i[0]) +#jvs
                #                             (self.n_BC[0] *  (conj(X)[1].dx(0)) +
                #                         self.n_BC[1] *   (conj(X)[1].dx(1))) * (fluc.u[1]-mean.u_forcing_i[1])) *self.ds(boundary_index)) #jvs

                # # Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
                # self.A_vf.add(-self.R*mean.nuTot*((self.n_BC[0] * (conj(X)[0].dx(0)) +
                #                         self.n_BC[1] * (conj(X)[0].dx(1) )) * (-mean.u_forcing_r[0]) +
                #                             (self.n_BC[0] * (conj(X)[1].dx(0) ) +
                #                         self.n_BC[1] * (conj(X)[1].dx(1))) * (-mean.u_forcing_r[1])) *self.ds(boundary_index))
    
    ## ---- Visc. 4: terms needed if we have three components to the mean velocity vector
    # These should be included in previous terms in tensor framework
    FLAG_TENS = True
    if param.Case.TransVelFluc:
        # -- > Previous implementation
        if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
            if FLAG_TENS:
                printWarning('--> Mom eq: Compressible + trans. vel not implemented in tensor framework. Currently relies on index notation.')
            
            # -- > Previous implementation
            self.A_vf.add(1j * self.R*Dx(mean.nuTot,i)*Dx(fluc.u[2],i)*conj(X)[2]*dx)       #III(1) (not 100% sure why only in case of dilatation -> needs checking)

            if not param.Case.m == 0:
                self.A_vf.add(1j * -self.m**2/self.R*mean.nuTot*fluc.u[2]*conj(X)[2]*dx)    #II(4)
                self.A_vf.add(1j * 2/3*self.m**2/self.R*mean.nuTot*fluc.u[2]*conj(X)[2]*dx) #IV(9)

                self.A_vf.add(-self.m*mean.nuTot*Dx(fluc.u[i],i)*conj(X)[2]*dx)			    #II(5)
                self.A_vf.add(-self.m/self.R*mean.nuTot*fluc.u[1]*conj(X)[2]*dx)			#II(6)
                self.A_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[i],i)*conj(X)[2]*dx)		#IV(7)
                self.A_vf.add(2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[1]*conj(X)[2]*dx)		#IV(8)

        if param.Case.AnalysisMode in ['Input-Output']:
            if FLAG_TENS:
                printWarning('--> Mom eq: Input-Output + trans. vel not implemented in tensor framework. Currently relies on index notation.')
                
            # -- > Previous implementation
            for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

                # Add physical boundary terms in imaginary part, where the forcing is applied...
                self.A_vf.add(1j * -self.R*mean.nuTot*(-conj(X) * (self.n_BC[0] *   (fluc.u[2].dx(0)) + self.n_BC[self.ThirdVelCompIndex] *   (fluc.u[2].dx(1)))) *self.ds(boundary_index))

                # Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
                self.A_vf.add(1j * -self.R * mean.nuTot * ((	self.n_BC[0] * (conj(X).dx(0)) +
                                                self.n_BC[self.ThirdVelCompIndex] * (conj(X).dx(1))) * (fluc.u[2] - mean.ut_forcing_i)) * self.ds(boundary_index))

            # Add stabilization terms on real part according to Baumann and Oden JFM 2016 vol 798
                self.A_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (conj(X).dx(0)) +
                                                self.n_BC[self.ThirdVelCompIndex] * (conj(X).dx(1))) * (-mean.ut_forcing_r)) * self.ds(boundary_index))
