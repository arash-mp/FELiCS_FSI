from ufl import (
				i,
				j,
				k,
				dx,
				Identity,
				)

def MomentumEqDilatationMax(self,fluc,X,mean,param):

	#lendim=param.mesh.geometry().dim()
	Iden=Identity(fluc.u.geometric_dimension())

	self.A_imag_vf.add( \
		-mean.nuTot*fluc.u[j].dx(i)*X[i].dx(j)*dx\
		+2./3. * mean.nuTot*fluc.u[k].dx(k)*Iden[i,j]*X[i].dx(j)*dx)
	self.A_imag_vf.add( \
		-fluc.nulam*mean.u[j].dx(i)*X[i].dx(j)*dx\
		+2./3. * fluc.nulam*mean.u[k].dx(k)*Iden[i,j]*X[i].dx(j)*dx)
