def addEnergyEq_LowMach(self,fluc,X,mean,param):
	from Equations.Energy_LowMach.EnergyEq_LowMachDivergenceU import EnergyEq_LowMachDivergenceU
	from Equations.Energy_LowMach.EnergyEq_LowMachLaplaceRho import EnergyEq_LowMachLaplaceRho
	from Equations.Energy_LowMach.EnergyEq_LowMachGradRhoGradRho import EnergyEq_LowMachGradRhoGradRho
	from FELiCS.functions import printDebug
	printDebug(param.debug,"Adding transport equation for density in all mesh internal directions")
	# rho^3 div( uFluc ) + 3 rho^2 rhoFluc div( uMean )
	EnergyEq_LowMachDivergenceU(self,fluc,X,mean,param)
	# 1/(Re Pr) * rhoMean Laplace( rhoFluc ) --> Integration by Parts is applied + 1/(Re Pr) * rhoFluc Laplace( rhoMean ) --> no IbP
	EnergyEq_LowMachLaplaceRho(self,fluc,X,mean,param)
	# 1/(Re Pr) * 4 * grad( rhoMean ) grad( rhoFluc )
	EnergyEq_LowMachGradRhoGradRho(self,fluc,X,mean,param)
	
	# the below function incorporates all of the above
#	from Equations.Energy_LowMach.tempAllTerms_MarioThesis import tempAllTerms_MarioThesis
#	tempAllTerms_MarioThesis(self,fluc,X,mean,param)
	

