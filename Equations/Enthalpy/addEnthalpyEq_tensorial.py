from ufl import (
                inner,
                dx,
                grad,
                div,
                conj,
                i
                )
from tensor_utils import *

def addEnthalpyEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    enthalpy conservation equation, in tensorial framework.
    
    NOTE:   The "FLAG_TENS" is used to activate the tensor
            framework in specific terms of the eq. It is possible
            to combine tensorial terms with index-notation ones
    '''

    from functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
    
    # Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED "mesh_dims" and "as_vector(x1, x2, 0.0)" to cartesian coords
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord_sys)  # HARDCODED FOR 2D perturbations so far
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)  # HARDCODED FOR 2D perturbations so far
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # CYL COORD IN FELiCS DEFINED AS (Z, R, PHI)!
        printWarning("Cylindrical coordinates for tensor species-cons not validated yet. Treat results with care.")
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord_sys)
        
    else:
        printError('Coord. syst not yet implemented in tensor framework.')
        
    # Definition of scalar terms in tensor framework
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys) 
    rho_mean_tens = Tensor(mean.rho, coord_sys)
    rho_fluc_tens = Tensor(fluc.rho, coord_sys)
    h_fluc_tens = Tensor(fluc.h, coord_sys)
    he_mean_tens = Tensor(mean.he, coord_sys)
    alpha_mean_tens = Tensor(mean.alpha, coord_sys)
    alpha_fluc_tens = Tensor(fluc.alpha, coord_sys)
    p_fluc_tens = Tensor(fluc.p, coord_sys)
    x_tens = Tensor(X, coord_sys)
    
    
    # ------------------------  Time derivative terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        self.B_vf.add((rho_mean_tens*h_fluc_tens*iConj(x_tens)).ufl_tens*coord_sys.J_hat*dx)
        self.B_vf.add((-1*p_fluc_tens*iConj(x_tens)).ufl_tens*coord_sys.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.B_vf.add( mean.rho*fluc.h*conj(X)*self.R*dx -fluc.p*conj(X)*self.R*dx)
    
    
    # ------------------------  Advection terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # Add volume integral of advection terms that remain after partial integration:
        self.A_vf.add((1j*iDiv(iConj(x_tens)*rho_mean_tens*u_mean_tens,-self.m)*h_fluc_tens).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add((1j*iDiv(iConj(x_tens)*rho_fluc_tens*u_mean_tens)*he_mean_tens).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add((1j*iDiv(iConj(x_tens)*rho_mean_tens*u_fluc_tens)*he_mean_tens).ufl_tens*coord_sys.J_hat*dx)
        
        # Add boundary integrals resulting from said partial integration:
        self.A_vf.add((-1j*iDot(rho_mean_tens*u_mean_tens*h_fluc_tens*iConj(x_tens),nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
        self.A_vf.add((-1j*iDot(rho_fluc_tens*u_mean_tens*he_mean_tens*iConj(x_tens),nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
        self.A_vf.add((-1j*iDot(rho_mean_tens*u_fluc_tens*he_mean_tens*iConj(x_tens),nbc_tens)).ufl_tens*coord_sys.J_hat*self.all_ds)
        
    else:
        # -- > Previous implementation
        #Add volume integral of advection terms that remain after partial integration:
        self.A_vf.add(1j * inner(div(conj(X)*mean.rho*mean.u*self.R),conj(fluc.h))*dx) # check whether conj(fluc.h) is right
        self.A_vf.add(1j * inner(div(conj(X)*mean.rho*fluc.u*self.R),mean.he)*dx)
        self.A_vf.add(1j * inner(div(conj(X)*fluc.rho*mean.u*self.R),mean.he)*dx)

        #Add boundary integrals resulting from said partial integration:
        self.A_vf.add(1j * -inner(conj(X)*mean.rho*mean.u*fluc.h*self.R,self.n_BC) * self.all_ds)
        self.A_vf.add(1j * -inner(conj(X)*mean.rho*fluc.u*mean.he*self.R,self.n_BC) * self.all_ds)
        self.A_vf.add(1j * -inner(conj(X)*fluc.rho*mean.u*mean.he*self.R,self.n_BC) * self.all_ds)

        #If coordinates are cylindrical add extra term for m not equal to zero:
        if param.Case.CoordinateSystem in ['Cylindrical'] and abs(self.m)>0:
            self.A_vf.add( X * self.m * mean.rho * mean.ut * fluc.h * dx)
        
    
        
        
    # ------------------------  Diffusion terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        self.A_vf.add(( -1j*iDot(alpha_mean_tens*iGrad(h_fluc_tens, self.m), iGrad(iConj(x_tens), -self.m)) ).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add(( -1j*iDot(alpha_fluc_tens*iGrad(he_mean_tens), iGrad(iConj(x_tens), -self.m)) ).ufl_tens*coord_sys.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.A_vf.add(1j * -inner(mean.alpha*grad(fluc.h),grad(X*self.R))*dx)	#Gleichung angepasst (CA)
        self.A_vf.add(1j * -fluc.alpha*mean.he.dx(i)*conj(X).dx(i)*dx)					#Alte Gleichung
        #self.A_imag_vf.add(-mean.alpha*fluc.viscoLaw*mean.he.dx(i)*X.dx(i)*dx)

        if param.Case.CoordinateSystem in ['Cylindrical']:
            self.A_vf.add(1j * mean.alpha*fluc.h.dx(1)*X*self.R*dx)
            if abs(self.m)>0:
                self.A_vf.add(1j * -mean.alpha*(self.m)*(self.m)/self.R*fluc.h*X*dx)

        #Debug test CA:
        #self.A_imag_vf.add(-inner(mean.alpha*mean.cp*grad(fluc.T),grad(X*self.R))*dx)	#Gleichung angepasst (CA)

        #if param.Case.CoordinateSystem in ['Cylindrical']:
        #	self.A_imag_vf.add(mean.alpha*mean.cp*fluc.T.dx(1)*X*self.R*dx)
        #	if abs(self.m)>0:
        #		self.A_imag_vf.add(-mean.alpha*mean.cp*(self.m)*(self.m)/self.R*fluc.T*X*dx)
        
    

