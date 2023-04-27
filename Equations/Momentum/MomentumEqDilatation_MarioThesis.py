from ufl import (
				i,
				j,
				dx,
				Dx,
				inner,
				grad,
				div,
				)
def MomentumEqDilatation(self,fluc,X,mean,param):
		# ************************************************************
		# The odl derivation neglected the transposed gradient tensor within the dissipative terms
		# The change within the linear equations also affects the dilatation terms
		# within Mario-Thesis this was changed, but it is not sophisticated
		DivergenceForm_real = div(fluc['u'])+fluc['u'][1]/self.R
		DivergenceForm_imag = self.m/self.R*fluc['ut']

		#x-comp
		self.A_imag_vf += 2*self.R/(3*mean['Re'])*X[0].dx(0)*DivergenceForm_real*dx
		if abs(self.m)>0:
			self.A_real_vf += -2*self.R/(3*mean['Re'])*X[0].dx(0)*DivergenceForm_imag*dx

		#r-comp
		self.A_imag_vf += 2*self.R/(3*mean['Re'])*X[1].dx(1)*DivergenceForm_real*dx
		if abs(self.m)>0:
			self.A_real_vf += -2*self.R/(3*mean['Re'])*X[1].dx(1)*DivergenceForm_imag*dx
