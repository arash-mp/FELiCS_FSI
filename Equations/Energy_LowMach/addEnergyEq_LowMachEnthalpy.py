from ufl import (
				as_vector,
				inner,
				derivative,
				grad,
				dx,
				i,
				)

def addEnergyEq_LowMachEnthalpy(self,fluc,X,mean,param, rup):
	#nonlinear space: rup
	ru  = as_vector((rup[0], rup[1]))
	rrho = rup[3]
###################
#derived variables

	print('Sutherland law for energy perturbation')
	M=28.949
	p=101300
	R=8314.4598/M

	Tm=p/R/mean['rho'] #mean T
	T=p/R/rrho #derived T

	Ts=170.672

	perturbT=(mean['T']+3*Ts)/(2*(mean['T']+Ts))*(-fluc['rho']/mean['rho']) #visco fluctuation factor(Sutherland)

	#derived sensible enthalpy
	h=mean['cp']*T+mean['he']-mean['cp']*Tm  #Explanation: he(variable)=Cp(MeanFlow)*T(variable)+he(MeanFlow)-Cp(MeanFlow)*T(MeanFlow)
	#h=mean['he']-mean['cp']*Tm  #Explanation: he(variable)=Cp(MeanFlow)*T(variable)+he(MeanFlow)-Cp(MeanFlow)*T(MeanFlow)
	#h=mean['cp']*T-mean['cp']*Tm  #Explanation: he(variable)=Cp(MeanFlow)*T(variable)+he(MeanFlow)-Cp(MeanFlow)*T(MeanFlow)
	#h=mean['he']-mean['cp']*Tm  #Explanation: he(variable)=Cp(MeanFlow)*T(variable)+he(MeanFlow)-Cp(MeanFlow)*T(MeanFlow)
	#after automatic derivative, we will get h'=Cp*T'

	#h0=h-(0.5*ru[i]*ru[i]) #used in diffusitivity term. Same in Max code,
	h0=h
###################
#non-linear form
	#F =-X*rrho*inner(ru,grad(h))*dx
	F =-X*mean['rho']*inner(ru,grad(h))*dx
	F += -X*rrho*inner(mean['u'],grad(mean['he']-mean['cp']*Tm))*dx
	#F =-X*rrho*inner(ru,grad(h))*dx
	#F =-X*rrho*inner(ru,grad(h))*dx
	#F = F- mean['alpha']*h0.dx(i)*X.dx(i)*dx
	##self.A_imag_vf += derivative(F, rup)

#time derivative

#	self.B_real_vf += self.rho*(mean['cp']*(Tm*(-fluc['rho']/mean['rho']))+inner(mean['u'],fluc['u']))*X*dx\
#	-fluc['p']*X*dx

# #add viscosity fluctuation
#	self.A_imag_vf += - perturbT*mean['alpha']*mean['he'].dx(i)*X.dx(i)*dx
