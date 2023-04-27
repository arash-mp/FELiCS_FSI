def addEnthalpyEq(self,fluc,X,mean,param):
	from Equations.Enthalpy.EnthalpyEqTimeDerivative import EnthalpyEqTimeDerivative
	from Equations.Enthalpy.EnthalpyEqAdvection import EnthalpyEqAdvection
	from Equations.Enthalpy.EnthalpyEqDiffusion import EnthalpyEqDiffusion
	from functions import printDebug
	printDebug(param.debug,"Adding transport equation for enthalpy in all mesh internal directions")
	EnthalpyEqTimeDerivative(self,fluc,X,mean)
	EnthalpyEqAdvection(self,fluc,X,mean,param)
	EnthalpyEqDiffusion(self,fluc,X,mean,param)
