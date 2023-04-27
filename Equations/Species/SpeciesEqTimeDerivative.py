from ufl import dx
def addSpeciesEqTimeDerivative(self,fluc,X,mean,species):
	self.B_real_vf.add( self.R*fluc.Y(species)*X*mean.rho*dx ) # Why is here no term fluc.rho * mean.Y
