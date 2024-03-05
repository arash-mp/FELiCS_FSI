from Equations.Species.SpeciesEqAdvection import addSpeciesEqAdvection
from Equations.Species.SpeciesEqTimeDerivative import addSpeciesEqTimeDerivative
from Equations.Species.SpeciesEqDiffusion import addSpeciesEqDiffusion
from Equations.Species.SpeciesEqSourceTerms import addSpeciesEqSourceTerms
from Equations.Species.SpeciesEqBCs import addSpeciesEqBCs
def addSpeciesEq(self,fluc,X,mean,species,param):
	addSpeciesEqAdvection(self,fluc,X,mean,species,param)
	addSpeciesEqTimeDerivative(self,fluc,X,mean,species)
	addSpeciesEqDiffusion(self,fluc,X,mean,species,param)
	addSpeciesEqSourceTerms(self,X,mean,species,param)
	addSpeciesEqBCs(self,fluc,X,mean,species,param)

