from fenics import as_vector,inner,derivative, grad, dx, Identity, div
from ufl import i,j,k

def addMomentumEq_LowMachEnthalpy(self,fluc,X,mean,param, rup):
	#nonlinear space: rup
	ru  = as_vector((rup[0], rup[1]))
	rp = rup[2]
	rrho = rup[3]



###################
#derived variables

	print('Sutherland law for perturbation')
	M=28.949
	p=101300
	R=8314.4598/M

	Tm=p/R/mean['rho'] #mean T
	T=p/R/rrho #derived T

	Ts=170.672

	perturbT=(Tm+3*Ts)/(2*(Tm+Ts))*(-fluc['rho']/mean['rho']) #visco fluctuation factor(Sutherland)



###################
#non-linear form
	lendim=param.mesh.geometry().dim()
	Iden=Identity(lendim)
	F   =  -( inner(grad(ru)*ru*rrho, X)*dx      \
	+ 1/mean['Re']*ru[i].dx(j)*X[i].dx(j)*dx + 1/mean['Re']*ru[j].dx(i)*X[i].dx(j)*dx\
	- 2/3 * 1/mean['Re']*ru[k].dx(k)*Iden[i,j]*X[i].dx(j)*dx-div(X)*rp*dx)
	self.A_imag_vf += derivative(F, rup)

#time derivative
	self.B_real_vf += self.R*(inner(fluc['u'],X))*self.rho *dx

# #add viscosity fluctuation
	self.A_imag_vf += \
		-perturbT* 1/mean['Re']*mean['u'][i].dx(j)*X[i].dx(j)*dx -perturbT*1/mean['Re']*mean['u'][j].dx(i)*X[i].dx(j)*dx\
		+perturbT* 2/3 * 1/mean['Re']*mean['u'][k].dx(k)*Iden[i,j]*X[i].dx(j)*dx