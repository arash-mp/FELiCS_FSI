# Third party libraries
import numpy as np


class viscosityHandler:
    """
    This class provides functions for calculating the viscosity.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution
    - meanFlowClass

    Private attributes:

    Protected attributes:

    Public attributes:

    """

    def __init__(self):
        pass

    def isMeanFlowClass(self):
        from meanFlowClass import meanFlowClass
        return isinstance(self, meanFlowClass)

    @property
    def MolViscModel(self):
        if self.isMeanFlowClass():
            return self._param.Case.MolViscModel
        else:
            return self._param.Case.MolViscPerturbModel

    @property
    def Viscosity(self):
        return self._param.Case.Mixture.Viscosity

    @property
    def As(self):
        return self.Viscosity['Constants']['As']

    @property
    def Ts(self):
        return self.Viscosity['Constants']['Ts']

    def getConstVisc(self):
        from dolfinx.fem import Function
        nulam = Function(self._FEMSpaces.P2)
        nulamVec = self._param.Case.MolVisc
        return nulam, nulamVec

    def getSutherlandViscMeanFlowClass(self, T, rho):
        from dolfinx.fem import Function, Constant
        muMol = self.As * (T.vector()[:] ** (3 / 2)) / (T.vector()[:] + self.Ts)
        nulam = Function(T.function_space())
        nulamVec = muMol / rho.vetor()[:]
        return nulam, nulamVec

    def getSutherlandViscFluctuationClass(self, mean, fieldDict, meanfieldDict,
                                          rho):
        from dolfinx.fem import Function, Constant
        if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != rho.shape[0]:
            nulam = self.As * fieldDict['T'] \
                    * (0.5
                       * meanfieldDict['T'].compute_vertex_values() ** (3 / 2)
                       + 1.5 * self.Ts
                       * meanfieldDict['T'].compute_vertex_values() ** 0.5) \
                    / ((meanfieldDict['T'].compute_vertex_values()
                        + self.Ts) ** 2)
        else:
            nulam = (self.As * fieldDict['T']
                     * (0.5 * mean.T ** (3 / 2) + 1.5 * self.Ts * mean.T ** 0.5)
                     / ((mean.T + self.Ts) ** 2))
        return nulam

    def getSutherlandMeanVisc(self, mean, meanfieldDict, rho):
        if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != rho.shape[0]:
            fluct = (meanfieldDict['T'].compute_vertex_values() + 3 * self.Ts) \
                    / (2 * (meanfieldDict['T'].compute_vertex_values()
                            + self.Ts)) \
                    * (-rho / meanfieldDict['rho'].compute_vertex_values())
            nulam = meanfieldDict['nulam'].compute_vertex_values() * fluct
        else:
            fluct = (mean.T + 3 * self.Ts) / (2 * (mean.T + self.Ts)) \
                    * (-rho / mean.rho)
            nulam = mean.nulam * fluct
        return nulam, fluct
