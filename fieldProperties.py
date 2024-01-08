# Third party libraries
from dolfinx.fem import (
    Constant,
    Function,
)
# Local Libraries and methods
from tensorUtils import (
                    iGrad,
                    iDiv,
                    iDot,
                    iIdentity,
                    iT,
                    Tensor,
                    )


class fieldProperties:
    """
    This class contains the property decorators for all field quantities by
    accessing the fieldDict dictionary that is used in fluctuationClass and
    meanFlowClass each.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution
    - meanFlowClass
    - meanFlowVertexValues

    Private attributes:

    Protected attributes:

    Public attributes:

    """

    def __init__(self):
        pass

    def isMeanFlowClass(self):
        from meanFlowClass import meanFlowClass
        return isinstance(self, meanFlowClass)

    def isMeanFlowVertexValuesClass(self):
        from meanFlowClass import meanFlowVertexValues
        return isinstance(self, meanFlowVertexValues)

    @property
    def alpha(self):
        if self.isMeanFlowClass:
            if 'alpha' in list(self._fieldDict.keys()):
                return  self._fieldDict['alpha']
                            
            else:
                return self._zeroField
        else:
            return Tensor(
                        self._fieldDict['alpha'],
                        self._coordinateSystem,
                        )

    @property
    def c(self):
        return self._fieldDict['c']

    @property
    def cp(self):
        return self._fieldDict['cp']

    def D(self, specie):
        if 'D_' + specie in list(self._fieldDict.keys()):
            return Tensor(
                        self._fieldDict['D_' + specie],
                        self._coordinateSystem,
                        )
        else:
            return Tensor(
                        self._zeroField,
                        self._coordinateSystem,
                        )
    @property
    def dQ(self):  # heat release
        return self._fieldDict['dQ']

    @property
    def fieldDict(self):
        if self.isMeanFlowClass:
            import copy
            Output = copy.copy(self._fieldDict)
            if '__hSpec' in dir(self):
                for key in list(self.__hSpec.keys()):
                    Output['hSpec_' + key] = self.__hSpec[key]
            return Output
        else:
            OutputDict = {}
            for key in list(self._fieldDict.keys()):
                if not isinstance(self._fieldDict[key], Function):
                    if key == 'u':
                        FEMSpace = self._FEMSpaces.FunctionSpaceVectorVelocity
                    else:
                        FEMSpace = self._FEMSpaces.P2
                    OutputDict[key] = project(self._fieldDict[key], FEMSpace)
                else:
                    OutputDict[key] = self._fieldDict[key]
            return OutputDict

    @fieldDict.setter
    def fieldDict(self, value):
        raise Exception('Properties of MeanFlow are not to be set after \
                        initialization!')

    @property
    def FieldNames(self):
        return self.__FieldNames

    @FieldNames.setter
    def FieldNames(self, value):
        raise Exception('Properties of MeanFlow are not to be set after \
                        initialization!')

    @property
    def fluc(self):
        return self._fluc

    def forcing_i(self, solution):
        return self._fieldDict[solution + '_forcing_i']

    def forcing_r(self, solution):
        return self._fieldDict[solution + '_forcing_r']

    @property
    def gamma(self):
        return self._fieldDict['gamma']

    @property
    def h(self):
        return self._fieldDict['h']

    @property
    def he(self):
        # To be generalized for all dimensions
        return Tensor(
                    self._fieldDict['he'],
                    self._coordinateSystem,
                    ) \
               + 0.5 * iDot(self.u, self.u)

    @property
    def hSpec(self):
        return self.__hSpec

    @property
    def meanflowFilename(self):
        return self._meanflowFilename

    @property
    def molarMass(self):
        return self._fieldDict['molarMass']

    @property
    def nulam(self):
        if 'nulam' in list(self._fieldDict.keys()):
            return self._fieldDict['nulam']
        else:
            return self._zeroField

    @property
    def nuTot(self):
        from dolfinx.fem import Function
        nuTot = Function(self._ScalarFunctionSpace)

        if 'nulam' in list(self._fieldDict.keys()):
            nuTot.vector[:] += self._fieldDict['nulam'].vector[:]
        if 'nuturb' in list(self._fieldDict.keys()):
            nuTot.vector[:] += self._fieldDict['nuturb'].vector[:]
        if 'nuSGS' in list(self._fieldDict.keys()):
            nuTot.vector[:] += self._fieldDict['nuSGS'].vector[:]
        return Tensor(
                    nuTot,
                    self._coordinateSystem,
                    )

    @property
    def p(self):
        if self.isMeanFlowClass():
            return Tensor(
                        self._fieldDict['p'],
                        self._coordinateSystem,
                        )
        else:
            if 'p' in self._transportedQuantities:
                return self._fieldDict['p']
            else:
                return Constant(0)

    @property
    def phi(self):
        if 'phi' in list(self._fieldDict.keys()):
            return self._fieldDict['phi']
        else:
            return self.__ZeroField

    @property
    def Q(self):  # heat release
        return self._fieldDict['Q']

    @property
    def reaction(self):
        return self.__reaction

    @property
    def rho(self):
        #if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
        if self.isMeanFlowClass():
            if 'rho' in list(self._fieldDict.keys()):
                return Tensor(
                            self._fieldDict['rho'],
                            self._coordinateSystem,
                            )
            else:
                from dolfinx.fem import Constant
                from petsc4py import PETSc
                mesh = self._fieldDict[list(self._fieldDict.keys())[0]].function_space.mesh
                return Constant(mesh, PETSc.ScalarType(1.0 + 0j))
        else:
            #if 'rho' in self._transportedQuantities:
            if 'rho' in list(self._fieldDict.keys()):
                return self._fieldDict['rho']
            else:
                return Tensor(
                            self._zeroField,
                            self._coordinateSystem,
                            )

    @property
    def rhou(self):
        return self._mean.rho * self.u + self.rho * self._mean.u

    def rhoY(self,species):
        return self._mean.rho * self.Y(species) + self.rho * self._mean.Y(species)

    @property
    def T(self):
        if self.isMeanFlowClass():
            from dolfinx.fem import Constant
            if 'T' in list(self._fieldDict.keys()):
                return Tensor(
                            self._fieldDict['T'],
                            self._coordinateSystem,
                            )
            else:
                return Constant(1) * self.__OneField
        else:
            return self._fieldDict['T']

    @property
    def tau(self):
        mean_nu = self._mean.nuTot
        mean_u = self._mean.u
        fluc_nu = self.nulam
        m = self._param.Case.m
        tau_out = mean_nu * iGrad(self.u,m) + \
                  fluc_nu * iGrad(mean_u)
        if not self._param.Case.SetOfEquations['Energy']['Equation'] == 'None':
            tau_out += iT(tau_out)
            tau_out += -2.0/3.0 * mean_nu * iDiv(self.u,m) * iIdentity(iGrad(self.u))
            tau_out += -2.0/3.0 * fluc_nu * iDiv(mean_u) * iIdentity(iGrad(self.u))
        return tau_out

    @property
    def Tb(self):
        if 'Tb' in list(self._fieldDict.keys()):
            return self._fieldDict['Tb']
        else:
            return Constant(1) * self.__OneField

    @property
    def Tu(self):
        if 'Tu' in list(self._fieldDict.keys()):
            return self._fieldDict['Tu']
        else:
            return Constant(1) * self.__OneField

    @property
    def Tm(self):
        return self._fieldDict['Tm']

    @property
    def responseDomain(self):
        return self._fieldDict['responseDomain']

    @property
    def forcingDomain(self):
        return self._fieldDict['forcingDomain']

    @property
    def u(self):
        if self.isMeanFlowClass():
            if 'u' in list(self._fieldDict.keys()):
                return Tensor(
                            self._fieldDict['u'],
                            self._coordinateSystem,
                            )
            else:
                return Tensor(
                            self._zeroVelocityField,
                            self._coordinateSystem,
                            )
        else:
            print(self._transportedQuantities)
            if 'u' in self._transportedQuantities:
                return  self._fieldDict['u']
            else:
                return Tensor(
                            self._zeroVelocityField,
                            self._coordinateSystem,
                            )

    @property
    def u_forcing_i(self):
        return self._fieldDict['u_forcing_i']

    @property
    def u_forcing_r(self):
        return self._fieldDict['u_forcing_r']

    @property
    def u_forcing(self):
        return Tensor(
                    self.u_forcing_r + 1j*self.u_forcing_i,
                    self._coordinateSystem,
                    )

    @property
    def ut(self):
        if self.isMeanFlowClass():
            if 'ut' in list(self._fieldDict.keys()):
                return self._fieldDict['u'][2]
            else:
                mesh = self._fieldDict[
                    list(self._fieldDict.keys())[0]].function_space.mesh
                return Constant(mesh, 0.0)
        else:
            if  self._param.Case.TransVelFluc:
                input('True')
                #return self._fieldDict['ut']
                return self._fieldDict['u'][2]
            else:
                input('False')
                return Constant(0)

    def Y(self, specie):
        if self.isMeanFlowClass():
            return Tensor(
                self._fieldDict[specie],
                self._coordinateSystem,
                )
        else:
            return self._fieldDict[specie]
                        
