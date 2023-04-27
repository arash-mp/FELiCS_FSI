

def getAlpha(param, mean, ClassDict, meanfieldDict, rho, isMeanFlowClass=True):
	import numpy as np
	from dolfinx.fem import Function,Constant
	
	if isMeanFlowClass:
		MolViscModel = param.Case.MolViscModel
	else:
		MolViscModel = param.Case.MolViscPerturbModel
	Pr 	= param.Case.Mixture.Pr
	Ts 	= param.Case.Mixture.Viscosity['Constants']['Ts']
	
	
	if MolViscModel == 	'Sutherland': 

		if isMeanFlowClass: 	# meanFlowClass
			alpha = Function(ClassDict['nulam'].function_space())
			alphaVector = ClassDict['nulam'].vector()[:]/Pr
			return alpha, alphaVector

		else: 		# FluctuationClass
			alpha = ClassDict['nulam'] / Pr
			return alpha	

	elif MolViscModel == 	'Sutherland mean':	# only used in FluctuationClass

		if 	(isinstance(mean.T, np.ndarray)
			 and mean.T.shape[0] != rho.shape[0]):	
			fluct =	(meanfieldDict['T'].compute_vertex_values() + 3*Ts) / (2 * (meanfieldDict['T'].compute_vertex_values() + Ts)) * (-rho / meanfieldDict['rho'].compute_vertex_values())
			alpha = meanfieldDict['nulam'].compute_vertex_values() * fluct
		else:
			fluct = (mean.T + 3*Ts) / (2 * (mean.T + Ts)) * (-rho / mean.rho)
			alpha = mean.alpha * fluct
		return alpha
