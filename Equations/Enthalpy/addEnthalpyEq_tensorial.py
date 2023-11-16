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
            
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(x)
    '''

    from functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
    
    # Coordinate system
    coord = self.coord_sys
    
    # Vector quantities
    if fluc.rhou.ufl_shape[0] == 2:
        u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord)
        u_m = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord)
        
    elif fluc.rhou.ufl_shape[0] == 3:
        u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord)
        u_m = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord)
    
    else:
        printError('--> NRJ-h eq: u has neither 2 or 3 dimensions: not implemented.') 
        
    # Scalar quantities
    rho_m = Tensor(mean.rho, coord)
    rho_f = Tensor(fluc.rho, coord)
    h_f = Tensor(fluc.h, coord)
    he_m = Tensor(mean.he, coord)
    alpha_m = Tensor(mean.alpha, coord)
    alpha_f = Tensor(fluc.alpha, coord)
    p_f = Tensor(fluc.p, coord)
    x_tens = Tensor(X, coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)
    
    
    # ------------------------  Time derivative terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        self.B_vf.add(( rho_m*h_f*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
        self.B_vf.add(( -1*p_f*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.B_vf.add( mean.rho*fluc.h*conj(X)*self.R*dx -fluc.p*conj(X)*self.R*dx)
    
    
    # ------------------------  Advection terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # Add volume integral of advection terms that remain after partial integration:
        self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_m*u_m,-self.m)*h_f ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_f*u_m)*he_m ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_m*u_f)*he_m ).ufl_tens*coord.J_hat*dx)
        # Add boundary integrals resulting from said partial integration:
        self.A_vf.add(( -1j*iDot(rho_m*u_m*h_f*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(rho_f*u_m*he_m*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(rho_m*u_f*he_m*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
        
    else:
        # -- > Previous implementation
        #Add volume integral of advection terms that remain after partial integration:
        self.A_vf.add(1j * inner(div(conj(X)*mean.rho*mean.u*self.R),conj(fluc.h))*dx) # check whether conj(fluc.h) is right
        self.A_vf.add(1j * inner(div(conj(X)*mean.rho*fluc.u*self.R),mean.he)*dx)
        self.A_vf.add(1j * inner(div(conj(X)*fluc.rho*mean.u*self.R),mean.he)*dx)
        # Add boundary integrals resulting from said partial integration:
        self.A_vf.add(1j * -inner(conj(X)*mean.rho*mean.u*fluc.h*self.R,self.n_BC) * self.all_ds)
        self.A_vf.add(1j * -inner(conj(X)*mean.rho*fluc.u*mean.he*self.R,self.n_BC) * self.all_ds)
        self.A_vf.add(1j * -inner(conj(X)*fluc.rho*mean.u*mean.he*self.R,self.n_BC) * self.all_ds)
        # If coordinates are cylindrical add extra term for m not equal to zero:
        if param.Case.CoordinateSystem in ['Cylindrical'] and abs(self.m)>0:
            self.A_vf.add( X * self.m * mean.rho * mean.ut * fluc.h * dx)
        
        
    # ------------------------  Diffusion terms
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        self.A_vf.add(( -1j*iDot(alpha_m*iGrad(h_f,self.m),iGrad(iConj(x_tens),-self.m)) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( -1j*iDot(alpha_f*iGrad(he_m),iGrad(iConj(x_tens),-self.m)) ).ufl_tens*coord.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.A_vf.add(1j * -inner(mean.alpha*grad(fluc.h),grad(X*self.R))*dx)	#Gleichung angepasst (CA)
        self.A_vf.add(1j * -fluc.alpha*mean.he.dx(i)*conj(X).dx(i)*dx)					#Alte Gleichung
        #self.A_imag_vf.add(-mean.alpha*fluc.viscoLaw*mean.he.dx(i)*X.dx(i)*dx)
        # If coordinates are cylindrical add extra term for m not equal to zero:
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
        
    

