from ufl import dx
from tensor_utils import *

def addEnergyPEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    energy conservation equation, in tensorial framework.
    The formulation is based on the total energy conservation
    for a compressible fluid, expressed in terms of pressure by
    substituting the Perfect gas law and continuity eq.
    The equation is implemented in PRIMITIVE variables. 
            
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(X)
    '''

    from functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding energy-p equation in all mesh internal directions")
    
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
    T_f = Tensor(fluc.T, coord)
    T_m = Tensor(mean.T, coord)
    p_f = Tensor(fluc.p, coord)
    p_m = Tensor(mean.p, coord)
    kappa_f = Tensor(fluc.kappa, coord)
    kappa_m = Tensor(mean.kappa, coord)
    gamma = mean.gamma
    x_tens = Tensor(X, coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord)
    
    
    # ------------------------  Time derivative terms
    # Volume term: -omega*p_f*conj(X)
    self.B_vf.add(( p_f*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
    
    
    # ------------------------  Advection terms
    # NOTE: "." denotes the dot product bellow
    # This includes two groups of terms:
    #   1)  -j*u.grad(p) = -j*(u_f.grad(p_m) + u_m.grad(p_f))
    #   2)  -j*gamma*p*div(u) =  -j*gamma*(p_f*div(u_m) + p_m*div(u_f))
    # both are integrated by parts in the weak form:
    # For 1)
    #   1.vol)  j*[div(conj(X)*u_f)*p_m + j*div(conj(X)*u_m)*p_f]*dx
    #   1.bc)   -j*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
    # Similarly for 2)
    #   2.vol)  j*gamma*[(grad(conj(X)*p_m)).u_f + (grad(conj(X)*p_f)).u_m]*dx
    #   2.bc)  -j*gamma*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
    # NOTE: the 1.bc) and 2.bc) terms are grouped together:
    #   bc)     -j*(1+gamma)*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
    
    # Add volume integral of pressure gradient term (1)
    self.A_vf.add(( 1j*iDiv(iConj(x_tens)*u_m,-self.m)*p_f ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(x_tens)*u_f)*p_m ).ufl_tens*coord.J_hat*dx)
    # Add volume integral of velocity divergence term (2)
    self.A_vf.add(( 1j*gamma*iDot(iGrad(iConj(X)*p_m,-self.m),u_f) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*gamma*iDot(iGrad(iConj(X)*p_f),u_m) ).ufl_tens*coord.J_hat*dx)
    # Add boundary integral of (1) and (2)
    self.A_vf.add(( -1j*(1+gamma)*iDot(u_m*p_f*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*(1+gamma)*iDot(u_f*p_m*iConj(x_tens),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)


    # ------------------------  Thermal diffusion term (Fourier law)
    # NOTE: kappa is the thermal conductivity
    # The term is: -j*(gamma-1)*div(kappa_m*grad(T_f) + kappa_f*grad(T_m))
    # which is integrated by parts,
    # Volume term:  j*(gamma-1)[(grad(conj(X)).(kappa_m*grad(T_f)) + (grad(conj(X)).(kappa_f*grad(T_m))]*dx
    self.A_vf.add(( 1j*(gamma-1)*iDot(iGrad(iConj(X),-self.m),kappa_m*iGrad(T_f,self.m)) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*(gamma-1)*iDot(iGrad(iConj(X),-self.m),kappa_f*iGrad(T_m)) ).ufl_tens*coord.J_hat*dx)
    # Boundary term:-j*(gamma-1)[(conj(X)*kappa_m*grad(T_f)).n_bc + (conj(X)*kappa_f*grad(T_m)).n_bc]*ds
    self.A_vf.add(( -1j*(gamma-1)*iDot(iConj(x_tens)*kappa_m*iGrad(T_f,self.m),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*(gamma-1)*iDot(iConj(x_tens)*kappa_f*iGrad(T_m),nbc_tens) ).ufl_tens*coord.J_hat*self.all_ds)
    
    
    # ------------------------  Viscous diffusion term
    # For now we neglect this term
        
    

