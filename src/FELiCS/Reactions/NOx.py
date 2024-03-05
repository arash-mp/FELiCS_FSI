

"""
"""


class NOx:
	def __init__(self, P):
		import pathlib
		from FELiCS.functions import loadCSV
		tablePath=str(pathlib.Path(__file__).parent.absolute())+'/NOxTable.csv'
		self.__Table=loadCSV(tablePath)



	def add_source_to_weak_form(self, weakform, dQ_threshold=None):
		"""
		Add a source term for the reaction and species to the weak form depending on the settings.
		dQ_threshold limits the reaction term to be applied only where the mean heat release rate
		is larger than dQ_threshold.
		"""
		from dolfinx.fem import (
								Function
								)
		from ufl import (
						dx,
		)
		from scipy import interpolate
		import numpy as np
		interpolationNO = interpolate.interp1d(self.__Table['phi'], self.__Table['omega_NO_pf'])
		interpolationNO2 = interpolate.interp1d(self.__Table['phi'], self.__Table['omega_NO2_pf'])
		omegaNO=Function(weakform.phi.function_space)
		print(np.max(weakform.phi.vector[:]))
		omegaNO.vector()[:]=interpolationNO(weakform.phi.vector[:])
		omegaNO2=Function(weakform.phi.function_space())
		omegaNO2.vector()[:]=interpolationNO2(weakform.phi.vector[:])
		print(len(omegaNO2.vector[:]))
		print(len(weakform.T.vector[:]))
		for i in range(len(omegaNO2.vector()[:])):
			if weakform.T.vector[i]<1000:
				omegaNO.vector[i]=0
				omegaNO2.vector[i]=0

		eq =  weakform.v_NO * omegaNO * dx
		eq += weakform.v_NO2 * omegaNO2 * dx
		return eq
