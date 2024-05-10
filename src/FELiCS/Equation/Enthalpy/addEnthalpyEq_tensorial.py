from ufl import dx
from FELiCS.Misc.tensorUtils import *

def addEnthalpyEq(self,fluc,X,mean,param):
    '''
    This function builds the weak form of the linearized
    enthalpy conservation equation, in tensorial framework.
    '''

    from FELiCS.Misc.functions import printDebug, printError, printWarning
    printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
    
    # Coordinate system
    J_hat = self._coordinateSystem.J_hat
    
    # ------------------------  Time derivative terms
    self.B_vf.add(( mean.rho*fluc.h*iConj(X) ).ufl_tens*J_hat*dx)
    self.B_vf.add(( -1*fluc.p*iConj(X) ).ufl_tens*J_hat*dx)
    

 
    # ------------------------  Advection terms
    # Add volume integral of advection terms that remain after partial integration:
    self.A_vf.add(( 1j*iDiv(iConj(X)*mean.rho*mean.u)*fluc.h ).ufl_tens*J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(X)*fluc.rho*mean.u)*mean.he ).ufl_tens*J_hat*dx)
    self.A_vf.add(( 1j*iDiv(iConj(X)*mean.rho*fluc.u)*mean.he ).ufl_tens*J_hat*dx)
    # Add boundary integrals resulting from said partial integration:
    self.A_vf.add(( -1j*iDot(mean.rho*mean.u*fluc.h*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(fluc.rho*mean.u*mean.he*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(mean.rho*fluc.u*mean.he*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
        
    # ------------------------  Diffusion terms
    self.A_vf.add(( -1j*iDot(mean.alpha*iGrad(fluc.h),iGrad(iConj(X))) ).ufl_tens*J_hat*dx)
    self.A_vf.add(( -1j*iDot(fluc.alpha*iGrad(mean.he),iGrad(iConj(X))) ).ufl_tens*J_hat*dx)

        
    

