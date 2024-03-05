from ufl import dx
def addSpeciesConservativeEqTimeDerivative(self,fluc,X,mean,species):
	self.B_real_vf.add( self.R*fluc.rhoY(species)*X*dx )
