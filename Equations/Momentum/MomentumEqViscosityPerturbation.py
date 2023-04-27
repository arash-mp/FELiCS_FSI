from ufl import (
				i,
				j,
				k,
				dx,
				Identity,
				)

def MomentumEqViscosityPerturbation(self,fluc,X,mean,param):

# ###################
# #derived variables
# 	print('Sutherland law for momentum perturbation')
# 	Ts=170.672
# 	perturbT=(mean.T+3*Ts)/(2*(mean.T+Ts))*(-fluc.rho/mean.rho) #visco fluctuation factor(Sutherland)

	# Get identity matrix with the size of the geometric dimension of the case
	Iden=Identity(fluc.rho.geometric_dimension())
	#lendim=param.mesh.geometry().dim()

# #add viscosity fluctuation
	self.A_imag_vf.add( \
		-fluc.viscoLaw* mean.nuTot*mean.u[i].dx(j)*X[i].dx(j)*dx\
		-fluc.viscoLaw* mean.nuTot*mean.u[j].dx(i)*X[i].dx(j)*dx\
		+fluc.viscoLaw* 2/3 * mean.nuTot*mean.u[k].dx(k)*Iden[i,j]*X[i].dx(j)*dx)
