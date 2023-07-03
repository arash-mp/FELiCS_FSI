from ufl import (
				i,
				j,
				k,
				Dx,
				dx,
				inner,
				outer,
				as_tensor,
				div,
				dot,
				grad,
				conj,
				Identity,
				)

import numpy as np

def MomentumEqAdvection(
	self,
	fluc,
	X,
	mean,
	param,
	):
	# Integrate convective terms in domain (integration by parts is applied)

	I=Identity( fluc.u.geometric_dimension() )
	self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * mean.u[j] ,k) * fluc.u[i] * I[k,j] * dx)
	self.A_vf.add(1j * Dx(self.R * conj(X[i]) * mean.rho * fluc.u[j] ,k) * mean.u[i] * I[k,j] * dx)
	self.A_vf.add(1j * Dx(self.R * conj(X[i]) * fluc.rho * mean.u[j] ,k) * mean.u[i] * I[k,j] * dx)

	# Integrate convective boundary terms (resulting from integration by parts)
	# Get convective flux F_conv
	F_conv_ij = ( mean.rho * mean.u[j] * fluc.u[i]\
		    + mean.rho * fluc.u[j] * mean.u[i]\
		    + fluc.rho * mean.u[j] * mean.u[i] )
	F_conv = as_tensor(F_conv_ij, (i, j))

	## Project flux normal to boundary (scalar product with n_BC) and add it to the momentum equation
	self.A_vf.add(1j * (- self.R * conj(X[i]) * mean.rho * mean.u[j] * fluc.u[i] * self.n_BC[k] * I[k,j]\
			    - self.R * conj(X[i]) * mean.rho * fluc.u[j] * mean.u[i] * self.n_BC[k] * I[k,j]\
			    - self.R * conj(X[i]) * fluc.rho * mean.u[j] * mean.u[i] * self.n_BC[k] * I[k,j])*self.all_ds)

	if param.Case.CoordinateSystem in ['Cylindrical']:
		self.A_vf.add(1j * - conj(X)[2] * mean.rho *mean.u[1] * fluc.u[2] * dx)
		self.A_vf.add(1j * - conj(X)[2] * mean.rho *fluc.u[1] * mean.u[2] * dx)
		self.A_vf.add(1j * - conj(X)[2] * fluc.rho *mean.u[1] * mean.u[2] * dx)

		self.A_vf.add(1j * conj(X)[1] * mean.u[2] * fluc.u[2] * mean.rho * dx\
					  + conj(X)[1] * fluc.u[2] * mean.u[2] * mean.rho * dx\
					  + conj(X)[1] * mean.u[2] * mean.u[2] * fluc.rho * dx)

		if not param.Case.m == 0:
			self.A_vf.add(mean.rho * conj(X)[i] * param.Case.m * mean.u[2] * fluc.u[i] * dx)

	if param.NumericalScheme in ['Discontinuous Galerkin']:
		#Lax-Friedrich-Flux coeffficient chosen as one
		LFFPrefactor=1
		# Calculate Lax-Friedrich fluxes
		C_LFF = LFFPrefactor * ( abs(dot(mean.u, self.n_BC)))
		F_u_ij = ( mean.u[i] * fluc.u[j] + mean.u[j] * fluc.u[i])
		F_u = as_tensor(F_u_ij, (i, j))  # Hier anstatt dieser und letzter Zeile F_conv von oben einfuegen?
		F_u_LFF = dot(avg(F_u), self.n_BC('+'))+ avg(C_LFF) * jump(fluc.u) / 2.0
		self.A_vf.add(1j * -dot(jump(X), F_u_LFF) * dS)
