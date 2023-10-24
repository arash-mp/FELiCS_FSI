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
    grad
)
from tensor_utils import (
    Tensor,
    CoordinateSystem,
    as_vector,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT
)
from functions import printWarning

def addMomentumEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    momentum conservation equation, in tensorial framework.
    The current version is based on the addMomentumEq.py function
    that was previously defined using the index notation of ufl.
    For now all terms from the momentum equation have been grouped 
    into a single file instead of being in serparate files.

    NOTE:
        (i) This version is hardcoded to CARTESIAN coordinates!
        (ii) It is an early attempt at integrating the tensor
            framework in FELiCS, so bugs/errors are expected.
        (iii) Support for discontinuous Galerkin has been dropped.
        (iv) Viscosity fluctuations are currently neglected.
        (v) This version is only for an INCOMPRESSIBLE flow, 
            and does not include density fluctuations.
    '''

    # ------------------------ Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED "mesh_dims" and "as_vector(x1, x2, 0.0)" to cartesian coords
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        x_tens = Tensor(as_vector((X[0], X[1], 0.0)), coord_sys)
        rho_mean_tens = Tensor(mean.rho, coord_sys)
        nutot_mean_tens = Tensor(mean.nuTot, coord_sys)
        nulam_fluc_tens = Tensor(fluc.nulam, coord_sys)
        rho_fluc_tens = Tensor(fluc.rho, coord_sys)
        p_fluc_tens = Tensor(fluc.p, coord_sys)
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)
        nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # HARDCODED FOR CYL COORD IN FELiCS, DEFINED AS (Z, R, PHI)!
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        x_tens = Tensor(as_vector((X[0], X[1], X[2])), coord_sys)
        rho_mean_tens = Tensor(mean.rho, coord_sys)
        nutot_mean_tens = Tensor(mean.nuTot, coord_sys)
        nulam_fluc_tens = Tensor(fluc.nulam, coord_sys)
        rho_fluc_tens = Tensor(fluc.rho, coord_sys)
        p_fluc_tens = Tensor(fluc.p, coord_sys)
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)
        nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)
        
    else:
        printWarning('Coord. syst not yet implemented in tensor framework.')

    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printWarning('Discontinuous Galerkin not implemented in tensorial framework.')

    if param.Case.CoordinateSystem =='Cylindrical':
        printWarning('Cyl. coord. not yet tested in tensor framework.')

    if not param.Case.m == 0:
        printWarning('m > 0 not yet implemented in tensor framework. \
                     For cyl. coord. it relies on previous ufl implementation.')

    if param.Case.AnalysisMode in ['Input-Output']:
        printWarning('Input-Outpout terms not yet implemented in tensor framework. \
                     It relies on previous ufl implementation.')


    # ------------------------ Time derivative term

    # -- > Previous implementation
    #self.B_vf.add(mean.rho*self.R*fluc.u[i]*conj(X)[i] *dx)
    # -- > Tensor implementation (not working for cyl coord)
    self.B_vf.add( (iDot(rho_mean_tens*u_fluc_tens, iConj(x_tens))).ufl_tens * coord_sys.J_hat * dx)


    # ------------------------ Convective terms

    ## Integrate convective terms in domain (integration by parts is applied)
    # -- > Previous implementation
    I=Identity( fluc.u.geometric_dimension() )
    #self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * mean.u[j] ,k) * fluc.u[i] * I[k,j] * dx)
    self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * fluc.u[j] ,k) * mean.u[i] * I[k,j] * dx)
    self.A_vf.add(1j * Dx(self.R * conj(X[i]) * fluc.rho * mean.u[j] ,k) * mean.u[i] * I[k,j] * dx)
    
    # -- > Tensor implementation which mimicks the previous implementation in FELiCS
    #self.A_vf.add( (1j * iDot(iDiv( iOuter(iConj(x_tens), rho_mean_tens*u_mean_tens) ), u_fluc_tens)).ufl_tens * coord_sys.J_hat * dx)
    #self.A_vf.add( (1j * iDot(iDiv( iOuter(iConj(x_tens), rho_mean_tens*u_fluc_tens) ), u_mean_tens)).ufl_tens * coord_sys.J_hat * dx)
    #self.A_vf.add( (1j * iDot(iDiv( iOuter(iConj(x_tens), rho_fluc_tens*u_mean_tens) ), u_mean_tens)).ufl_tens * coord_sys.J_hat * dx)
    
    # -- > Tensor implementation derived by hand (deosn't seem to work)
    self.A_vf.add( (1j * iInner(iGrad( iConj(x_tens) ), rho_mean_tens*iOuter( u_fluc_tens, u_mean_tens )) ).ufl_tens * coord_sys.J_hat * dx)
    #self.A_vf.add( (1j * iInner(iGrad( iConj(x_tens) ), rho_mean_tens*iOuter( u_mean_tens, u_fluc_tens )) ).ufl_tens * coord_sys.J_hat * dx)
    #self.A_vf.add( (1j * iInner(iGrad( iConj(x_tens) ), rho_fluc_tens*iOuter( u_mean_tens, u_mean_tens )) ).ufl_tens * coord_sys.J_hat * dx)
    
    # NOTE: The last term above has not been validated yet.
    
    if param.Case.CoordinateSystem in ['Cylindrical']:
        # The terms below should no longer be necessary in tensor framework
        #self.A_vf.add(1j * - conj(X)[2] * mean.rho *mean.u[1] * fluc.u[2] * dx)
        #self.A_vf.add(1j * - conj(X)[2] * mean.rho *fluc.u[1] * mean.u[2] * dx)
        #self.A_vf.add(1j * - conj(X)[2] * fluc.rho *mean.u[1] * mean.u[2] * dx)

        #self.A_vf.add(1j * conj(X)[1] * mean.u[2] * fluc.u[2] * mean.rho * dx)
        #self.A_vf.add(1j * conj(X)[1] * fluc.u[2] * mean.u[2] * mean.rho * dx)
        #self.A_vf.add(1j * conj(X)[1] * mean.u[2] * mean.u[2] * fluc.rho * dx)

        if not param.Case.m == 0:
            self.A_vf.add(mean.rho * conj(X)[i] * param.Case.m * mean.u[2] * fluc.u[i] * dx)

    ## Project flux normal to boundary (scalar product with n_BC)
    # and add it to the momentum equation
    # -- > Previous implementation
    self.A_vf.add(1j * (- self.R * conj(X[i]) * mean.rho * mean.u[j] * fluc.u[i] * self.n_BC[k] * I[k,j])*self.all_ds)
    self.A_vf.add(1j * (- self.R * conj(X[i]) * mean.rho * fluc.u[j] * mean.u[i] * self.n_BC[k] * I[k,j])*self.all_ds)
    self.A_vf.add(1j * (- self.R * conj(X[i]) * fluc.rho * mean.u[j] * mean.u[i] * self.n_BC[k] * I[k,j])*self.all_ds)
    
    # -- > Tensor implementation which mimicks the previous implementation in FELiCS
    #self.A_vf.add( (-1j * rho_mean_tens * iDot(u_mean_tens, nbc_tens) * iDot(u_fluc_tens, iConj(x_tens)) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * rho_mean_tens*iDot(u_fluc_tens, nbc_tens) * iDot(u_mean_tens, iConj(x_tens)) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * rho_fluc_tens*iDot(u_mean_tens, nbc_tens) * iDot(u_mean_tens, iConj(x_tens)) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    
    # -- > Tensor implementation derived by hand
    #self.A_vf.add( (-1j * iDot( iDot( (rho_mean_tens*iOuter(u_fluc_tens, u_mean_tens)), iConj(x_tens) ), nbc_tens) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * iDot( iDot( (rho_mean_tens*iOuter(u_mean_tens, u_fluc_tens)), iConj(x_tens) ), nbc_tens) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * iDot( iDot(rho_fluc_tens * iOuter(u_mean_tens, u_mean_tens), iConj(x_tens)), nbc_tens) ).ufl_tens*coord_sys.J_hat*self.all_ds)
    
    #self.A_vf.add( (-1j * iDot(iDot(iConj(x_tens), rho_mean_tens * iOuter(u_fluc_tens, u_mean_tens)), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * iDot(iDot(iConj(x_tens), rho_mean_tens * iOuter(u_mean_tens, u_fluc_tens)), nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
    #self.A_vf.add( (-1j * iDot(iDot(iConj(x_tens), rho_fluc_tens * iOuter(u_mean_tens, u_mean_tens)), nbc_tens) ).ufl_tens*coord_sys.J_hat*self.all_ds)


    # ------------------------ Pressure gradient terms

    int_by_parts = True
    if int_by_parts:
        # Integrate pressure gradient boundary terms (resulting from integration by parts)
        # -- > Previous implementation
        #self.A_vf.add(1j * Dx(conj(X[i])*self.R,j)*fluc.p * I[i,j] * dx)
        #self.A_vf.add(1j * -self.R*fluc.p*self.n_BC[j]*conj(X[i])* I[i,j] * self.all_ds)
        # -- > Tensor implementation (not working for cyl coord)
        self.A_vf.add( (1j * p_fluc_tens*iDiv(iConj(x_tens)) ).ufl_tens * coord_sys.J_hat * dx)
        self.A_vf.add( (-1j * iDot(p_fluc_tens*iConj(x_tens), nbc_tens))\
                    .ufl_tens * coord_sys.J_hat * self.all_ds)

    else:
        printWarning('Pressure term without integration by parts not \
                     yet validated in tensor framework.')
        # -- > Previous implementation
        #self.A_vf.add(1j * -self.R*inner(grad(fluc.p), X)*dx)
        # -- > Tensor implementation (not working for cyl coord)
        self.A_vf.add( (-1j * iDot(iGrad(p_fluc_tens), iConj(x_tens)) )\
                      .ufl_tens * coord_sys.J_hat * dx)

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

    # ---- Visc. 1: viscous terms for incompressible flow
    # -- > Previous implementation
    #self.A_vf.add(1j * -self.R*mean.nuTot*Dx(fluc.u[i],j)*Dx(conj(X)[i],j)*dx )
    #self.A_vf.add(1j * -fluc.nulam*mean.u[i].dx(j)*conj(X)[i].dx(j)*dx )
    # -- > Tensor implementation (not working for cyl coord)
    self.A_vf.add( (-1j * nutot_mean_tens * iInner(iGrad(u_fluc_tens), iGrad(iConj(x_tens))) )\
                .ufl_tens * coord_sys.J_hat * dx)
    self.A_vf.add( (-1j * nulam_fluc_tens * iInner(iGrad(u_mean_tens), iGrad(iConj(x_tens))) )\
                .ufl_tens * coord_sys.J_hat * dx)
    # NOTE: The term with nulam_fluc has not been validated yet

    # ---- Visc. 3: additional viscous terms for incompressible flow in cyl. coord
    #if param.Case.CoordinateSystem =='Cylindrical':									#is this necessary? (CA)
        #self.A_vf.add(1j * -mean.nuTot/self.R*conj(X)[2]*fluc.u[2]*dx)#
        #self.A_vf.add(1j * -mean.nuTot.dx(1)*conj(X)[2]*fluc.u[2]*dx)
        #self.A_vf.add(1j * -mean.nuTot*conj(X)[1]/self.R*fluc.u[1]*dx)#
    if not param.Case.m == 0:
        self.A_vf.add(-2*mean.nuTot/self.R*self.m*conj(X)[2]*fluc.u[1]*dx)#
        self.A_vf.add(-conj(X)[2]*mean.nuTot.dx(0)*self.m*fluc.u[0]*dx)
        #self.A_vf.add(1j * -X*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)							#III  -> A_real ??? (CA)
        self.A_vf.add(-conj(X)[2]*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)
        self.A_vf.add(1j * -mean.nuTot*self.m**2/self.R*inner(fluc.u, X)*dx)#
        #self.A_vf.add(1j * -mean.nuTot*self.m**2/self.R*conj(X)[i]*fluc.u[i]*dx)# this line replaces the above which uses inner
        self.A_vf.add(2*mean.nuTot*self.m/self.R*conj(X)[1]*fluc.u[2]*dx)#
    
    # ---- Visc. 5: viscous BC terms
    if param.Case.AnalysisMode in ['Input-Output']:
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
            self.A_vf.add(1j * -self.R*mean.nuTot*((self.n_BC[0] * (conj(X)[0].dx(0)) +
                                       self.n_BC[1] *   (conj(X)[0].dx(1))) * (fluc.u[0]-mean.u_forcing_i[0]) +#jvs
                                         (self.n_BC[0] *  (conj(X)[1].dx(0)) +
                                       self.n_BC[1] *   (conj(X)[1].dx(1))) * (fluc.u[1]-mean.u_forcing_i[1])) *self.ds(boundary_index)) #jvs

            # Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
            self.A_vf.add(-self.R*mean.nuTot*((self.n_BC[0] * (conj(X)[0].dx(0)) +
                                       self.n_BC[1] * (conj(X)[0].dx(1) )) * (-mean.u_forcing_r[0]) +
                                         (self.n_BC[0] * (conj(X)[1].dx(0) ) +
                                       self.n_BC[1] * (conj(X)[1].dx(1))) * (-mean.u_forcing_r[1])) *self.ds(boundary_index))
    
    # ---- Visc. 6: other stuff
    if  param.Case.TransVelFluc:
        if not param.Case.SetOfEquations['Energy']['Equation'] == 'None':
            self.A_vf.add(1j * self.R*Dx(mean.nuTot,i)*Dx(fluc.u[2],i)*conj(X)[2]*dx)						#III(1) (not 100% sure why only in case of dilatation -> needs checking)

            if not param.Case.m == 0:
                self.A_vf.add(1j * -self.m**2/self.R*mean.nuTot*fluc.u[2]*conj(X)[2]*dx)					#II(4)
                self.A_vf.add(1j * 2/3*self.m**2/self.R*mean.nuTot*fluc.u[2]*conj(X)[2]*dx)					#IV(9)

                self.A_vf.add(-self.m*mean.nuTot*Dx(fluc.u[i],i)*conj(X)[2]*dx)						#II(5)
                self.A_vf.add(-self.m/self.R*mean.nuTot*fluc.u[1]*conj(X)[2]*dx)						#II(6)
                self.A_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[i],i)*conj(X)[2]*dx)				#IV(7)
                self.A_vf.add(2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[1]*conj(X)[2]*dx)				#IV(8)




        if param.Case.AnalysisMode in ['Input-Output']:
            for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

                # Add physical boundary terms in imaginary part, where the forcing is applied...
                self.A_vf.add(1j * -self.R*mean.nuTot*(-conj(X) * (self.n_BC[0] *   (fluc.u[2].dx(0)) + self.n_BC[self.ThirdVelCompIndex] *   (fluc.u[2].dx(1)))) *self.ds(boundary_index))

                # Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
                self.A_vf.add(1j * -self.R * mean.nuTot * ((	self.n_BC[0] * (conj(X).dx(0)) +
                                                self.n_BC[self.ThirdVelCompIndex] * (conj(X).dx(1))) * (fluc.u[2] - mean.ut_forcing_i)) * self.ds(boundary_index))

            # Add stabilization terms on real part according to Baumann and Oden JFM 2016 vol 798
                self.A_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (conj(X).dx(0)) +
                                                self.n_BC[self.ThirdVelCompIndex] * (conj(X).dx(1))) * (-mean.ut_forcing_r)) * self.ds(boundary_index))
