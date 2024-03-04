from ufl import (
				i,
				j,
				dx,
				inner,
				as_tensor,
				div,
				dot,
				grad,
				Dx,
				Identity,
				conj,
				)

def MomentumEqPressureGradient(self,fluc,X,mean,param):
	IntegrationByParts=True
	# Integrate pressure gradient in domain
	print(fluc.u)
	print(fluc.p)
	print(X)
	I=Identity( fluc.u.geometric_dimension() )
	if IntegrationByParts:
		#self.A_vf.add(1j * inner(div(X*self.R),fluc.p)*dx)
		self.A_vf.add(1j * Dx(conj(X[i])*self.R,j)*fluc.p * I[i,j] * dx)
		# Integrate pressure gradient boundary terms (resulting from integration by parts)
		self.A_vf.add(1j * -self.R*fluc.p*self.n_BC[j]*conj(X[i])* I[i,j] * self.all_ds)
	else:
		self.A_vf.add(1j * -self.R*inner(grad(fluc.p), X)*dx)

	if not param.Case.m == 0:
		self.A_vf.add(+conj(X)[2]*self.m*fluc.p*dx)


	if param.NumericalScheme in ['Discontinuous Galerkin']:
		if param.Case.CoordinateSystem in ['Cylindrical']:
			printError("Discontinuous Elements not implemented for cylindrical coordinates")
		#Lax-Friedrich-Flux coeffficient chosen as one
		LFFPrefactor=1
		# Calculate Lax-Friedrich fluxes
		C_LFF = LFFPrefactor * ( abs(dot(mean.u, self.n_BC)))
		Iden = Identity(2)
		F_p_ij = Iden[i,j] * fluc.p
		F_p = as_tensor(F_p_ij, (i, j))
		F_p_LFF = dot(avg(F_p), self.n_BC('+'))+ avg(C_LFF) * jump(fluc.u) / 2.0
		self.A_vf.add(1j * -aot(jump(X), F_p_LFF) * dS)