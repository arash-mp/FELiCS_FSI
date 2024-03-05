from Equations.Mass.MassEqAdvection import MassEqAdvection
def addMassEq(self,fluc,X,mean,param):
	MassEqAdvection(self,fluc,X,mean,param)
	if 'rho' in param.Case.getTransportedQuantityList():
		from Equations.Mass.MassEqTimeDerivative import MassEqTimeDerivative
		MassEqTimeDerivative(self,fluc,X,mean,param)
