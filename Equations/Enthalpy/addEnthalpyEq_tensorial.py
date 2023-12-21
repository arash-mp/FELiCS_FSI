from ufl import dx
from tensorUtils import *

def addEnthalpyEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    enthalpy conservation equation, in tensorial framework.
            
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(X)
    '''

    from functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
    
    # Coordinate system
    coord = self.coord_sys
    
    # Vector quantities
    if fluc.rhou.ufl_shape[0] == 2:
        printDebug(param.debug, '--> Enthalpy cons. eq: fluctuations and mean flow 2D')
        u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord)
        u_m = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord)
        
    elif fluc.rhou.ufl_shape[0] == 3:
        printDebug(param.debug, '--> Enthalpy cons. eq: fluctuations and mean flow 3D')
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
    x_tens = Tensor(
                    X,
                    coord,
                    containsTestFunction = True,
                    )
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)
    
    
    # ------------------------  Time derivative terms
    self.B_vf.add(( rho_m*h_f*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
    self.B_vf.add(( -1*p_f*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
    
    
    # ------------------------  Advection terms
    # Add volume integral of advection terms that remain after partial integration:
    self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_m*u_m,self.m)*h_f ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_f*u_m)*he_m ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(x_tens)*rho_m*u_f)*he_m ).ufl_tens*coord.J_hat*dx)
    # Add boundary integrals resulting from said partial integration:
    self.A_vf.add(( -1j*iDot(rho_m*u_m*h_f*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(rho_f*u_m*he_m*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(rho_m*u_f*he_m*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
        
        
    # ------------------------  Diffusion terms
    self.A_vf.add(( -1j*iDot(alpha_m*iGrad(h_f,self.m),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( -1j*iDot(alpha_f*iGrad(he_m),iGrad(iConj(x_tens),self.m)) ).ufl_tens*coord.J_hat*dx)

        
    

