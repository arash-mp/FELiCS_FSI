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

    from Equations.Enthalpy.EnthalpyEqAdvection import EnthalpyEqAdvection
    from Equations.Enthalpy.EnthalpyEqDiffusion import EnthalpyEqDiffusion
    from functions import printDebug
    printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
    
    
    # ------------------------  Time derivative terms
    self.B_vf.add( mean.rho*fluc.h*conj(X)*self.R*dx-fluc.p*conj(X)*self.R*dx)
    
    
    # ------------------------  Advection terms
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

