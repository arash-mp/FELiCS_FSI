from fenics import as_vector,inner,derivative, grad, dx, div,dot
from ufl import i,j,k

def addMassEq_LowMachEnthalpy(self,fluc,X,mean,param, rup):
	#nonlinear space: rup
	ru  = as_vector((rup[0], rup[1]))
	rrho = rup[3]

###################
#non-linear form
	F   =  -(dot(grad(rrho),ru)*X*dx \
		+ X*rrho*div(ru)*dx)
	self.A_imag_vf += derivative(F, rup)



#time derivative
	self.B_real_vf += fluc['rho']*self.R*X*dx

