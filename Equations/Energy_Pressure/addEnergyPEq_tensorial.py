from ufl import dx
from tensorUtils import *

def addEnergyPEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    energy conservation equation, in tensorial framework.
    The formulation is based on the total energy conservation
    for a compressible fluid, expressed in terms of pressure by
    substituting the Perfect gas law and continuity eq.
    The equation is implemented in PRIMITIVE variables. 
    '''

    from functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding energy-p equation in all mesh internal directions")
    
    # ------------------------ Define the tensorial operators
    J_hat = self._coordinateSystem.J_hat
    
    # ------------------------  Time derivative terms
    # Volume term: -omega*p_f*conj(X)
    self.B_vf.add(( fluc.p*iConj(X) ).ufl_tens*J_hat*dx)
    
    
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
    self.A_vf.add(( 1j*iDiv(iConj(X)*mean.u)*fluc.p ).ufl_tens*J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(X)*fluc.u)*mean.p ).ufl_tens*J_hat*dx)
    # Add volume integral of velocity divergence term (2)
    self.A_vf.add(( 1j*mean.gamma*iDot(iGrad(iConj(X)*mean.p),fluc.u) ).ufl_tens*J_hat*dx)
    self.A_vf.add(( 1j*mean.gamma*iDot(iGrad(iConj(X)*fluc.p),mean.u) ).ufl_tens*J_hat*dx)
    # Add boundary integral of (1) and (2)
    self.A_vf.add(( -1j*(1+mean.gamma)*iDot(mean.u*fluc.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
    self.A_vf.add(( -1j*(1+mean.gamma)*iDot(fluc.u*mean.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)


    # ------------------------  Thermal diffusion term (Fourier law)
    # NOTE: kappa is the thermal conductivity
    # The term is: -j*(gamma-1)*div(kappa_m*grad(T_f) + kappa_f*grad(T_m))
    # which is integrated by parts,
    # TEMPORARY: use constant Pr definition (This should move to fieldProperty or a handler)
    kappa_m = mean.nuTot*mean.cp/mean.Pr
    kappa_f = fluc.nuTot*mean.cp/mean.Pr
    # Volume term:  
    #   j*(gamma-1)[(grad(conj(X)).(kappa_m*grad(T_f)) + (grad(conj(X)).(kappa_f*grad(T_m))]*dx
    self.A_vf.add(( 1j*(mean.gamma-1)*iDot(iGrad(iConj(X)),kappa_m*iGrad(fluc.T)) ).ufl_tens*J_hat*dx)
    self.A_vf.add(( 1j*(mean.gamma-1)*iDot(iGrad(iConj(X)),kappa_f*iGrad(mean.T)) ).ufl_tens*J_hat*dx)
    # Boundary term:
    #   -j*(gamma-1)[(conj(X)*kappa_m*grad(T_f)).n_bc + (conj(X)*kappa_f*grad(T_m)).n_bc]*ds
    self.A_vf.add(( -1j*(mean.gamma-1)*iDot(iConj(X)*kappa_m*iGrad(fluc.T),self.n) ).ufl_tens*J_hat*self.all_ds)
    self.A_vf.add(( -1j*(mean.gamma-1)*iDot(iConj(X)*kappa_f*iGrad(mean.T),self.n) ).ufl_tens*J_hat*self.all_ds)
    
    
    # ------------------------  Viscous diffusion term
    # For now we neglect this term
        
    

