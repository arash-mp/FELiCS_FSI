from ufl import (
				dx,
				conj,
				)
def addSpeciesEqTimeDerivative(self,fluc,X,mean,species):
	self.B_vf.add( self.R*fluc.Y(species)*conj(X)*mean.rho*dx )

