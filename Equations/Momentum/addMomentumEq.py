from ufl import(
				conj,
)
def addMomentumEq(self,fluc,X,mean,param):
	from Equations.Momentum.MomentumEqTimeDerivative import MomentumEqTimeDerivative
	from Equations.Momentum.MomentumEqAdvection import MomentumEqAdvection
	from Equations.Momentum.MomentumEqPressureGradient import MomentumEqPressureGradient
	from Equations.Momentum.MomentumEqDiffusionMario2 import MomentumEqDiffusion
	from Equations.Momentum.MomentumEqDilatationMario2 import MomentumEqDilatation
	from functions import printDebug
	printDebug(param.debug,"Adding transport equation for momentum in all mesh internal directions")
	MomentumEqTimeDerivative(self,fluc,conj(X),mean)
	MomentumEqAdvection(self,fluc,X,mean,param)
	MomentumEqPressureGradient(self,fluc,X,mean,param)
	MomentumEqDiffusion(self,fluc,conj(X),mean,param)
	# if needed, add dilatation to momentum equation
	#if param.FlowMode in ['LowMach', 'Reacting' ]:
	#	MomentumEqDilatation(self,fluc,X,mean,param)

	##Add viscosity perturbation
	#from Equations.Momentum.MomentumEqViscosityPerturbation import MomentumEqViscosityPerturbation
	#MomentumEqViscosityPerturbation(self,fluc,X,mean,param)
	#if param.Case.HeatTransfer or param.Case.Compressible:
	#	#Add dilatation in Max's way
	#	from Equations.Momentum.MomentumEqDilatationMax import MomentumEqDilatationMax
	#	MomentumEqDilatationMax(self,fluc,X,mean,param)

