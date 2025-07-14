# Third party libraries
from dolfinx.fem import (
    Constant,
    Function
)
# Local Libraries and methods
from FELiCS.Misc.tensorUtils import (
                    iGrad,
                    iDiv,
                    iDot,
                    iIdentity,
                    iT,
                    Tensor,
                    )
from FELiCS.Fields.Field import Field

from FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

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

    def isMeanFlowClass(self):
        from FELiCS.Fields.meanFlowClass import meanFlowClass
        return isinstance(self, meanFlowClass)

    def isMeanFlowVertexValuesClass(self):
        from FELiCS.Fields.meanFlowClass import meanFlowVertexValues
        return isinstance(self, meanFlowVertexValues)

    @property
    def alpha(self):
        if self.isMeanFlowClass():
            if 'alpha' in list(self._fieldDict.keys()):
                return  self._fieldDict['alpha'].getTensor()
            else:
                return self._zeroField.getTensor()
        if self.isMeanFlowVertexValuesClass():
            if 'alpha' in list(self._fieldDict.keys()):
                return self._fieldDict['alpha']
            else:
                return self._zeroField
        else:
            return self._fieldDict['alpha']

    @property
    def c(self):
        return self._fieldDict['c']

    @property
    def cp(self):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            return self._fieldDict['cp'].getTensor()
        else:       
            return self._fieldDict['cp']

    def D(self, specie):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            if 'D_' + specie in list(self._fieldDict.keys()):
                return self._fieldDict['D_' + specie].getTensor()
            else:
                return self._zeroField.getTensor()
        else:
            if 'D_' + specie in list(self._fieldDict.keys()):
                return self._fieldDict['D_' + specie]
            else:
                return self._zeroField


    @property
    def dQ(self):  # heat release
        return self._fieldDict['dQ']

    @property
    def fieldDict(self):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
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
                    OutputDict[key] = project(self._fieldDict[key], FEMSpace) # NOTE: (Simon) not sure what this whould be
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

    def forcing(self,solution):
        result = self.forcing_r(solution)+self.forcing_i(solution)*1j
        return result.getTensor()

    @property
    def forcingDomain(self):
        return self._fieldDict['forcingDomain']

    @property
    def gamma(self):    # Heat capacity ratio
        return self._fieldDict['gamma'].getTensor()
        
    @property
    def Pr(self):    # Prandtl number 
        return self._fieldDict['Pr'].getTensor()
    
    @property
    def UnitT(self):    # Real unit number in tensor form
        from dolfinx.fem import Constant
        from petsc4py import PETSc
        mesh = self._fieldDict[list(self._fieldDict.keys())[0]].space.mesh
        return Tensor(
            Constant(mesh, PETSc.ScalarType(1.0 + 0j)),
            self._coordinateSystem,
            )

    @property
    def h(self):
        return self._fieldDict['h']

    @property
    def he(self):
        return self._fieldDict['he'].getTensor() \
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
        if self.isMeanFlowClass():
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam'].getTensor()
            else:
                return self._zeroField.getTensor()
        elif self.isMeanFlowVertexValuesClass():
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam']
            else:
                return self._zeroField
        
        else:
            if 'nulam' in list(self._fieldDict.keys()):
                return self._fieldDict['nulam']
            else:
                return self._zeroField

    @property
    def nuTot(self):
        from dolfinx.fem import Function
        nuTot = Field(self._ScalarFunctionSpace, self._mesh)

        if 'nulam' in list(self._fieldDict.keys()):
            nuTot += self._fieldDict['nulam']
        if 'nuturb' in list(self._fieldDict.keys()):
            nuTot += self._fieldDict['nuturb']
        if 'nuSGS' in list(self._fieldDict.keys()):
            nuTot += self._fieldDict['nuSGS']
        return nuTot.getTensor()

    @property
    def p(self):
        if self.isMeanFlowClass():
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p'].getTensor()
            else:
                return self._zeroField.getTensor()
            
        elif self.isMeanFlowVertexValuesClass():
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p']
            else:
                return self._oneField

        else:
            if 'p' in list(self._fieldDict.keys()):
                return self._fieldDict['p']
            else:
                return self._zeroField
            
    @property
    def phi(self):
        if 'phi' in list(self._fieldDict.keys()):
            return self._fieldDict['phi']
        else:
            return self._zeroField

    @property
    def RR_prefactor(self):
        if self.isMeanFlowClass():
            if 'RR_prefactor' in list(self._fieldDict.keys()):
                return self._fieldDict['RR_prefactor'].getTensor()
        else:
            return self._fieldDict['RR_prefactor']

    @property
    def Q(self):  # heat release
        return self._fieldDict['Q']
    
    @property
    def R_spe(self):  # Specific gas constant
        return self._fieldDict['R_spe'].getTensor()

    @property
    def reaction(self):
        return self.__reaction
    
    @property
    def spg(self):  # Sponge region term
        return self._fieldDict['spg'].getTensor()

    @property
    def rho(self):
        if self.isMeanFlowClass():
            # Comment from Sophie: I added rho to the quantities to read in as default, s.t. a variable density
            # without rho as fluctuation variable is possible ("cold flow"). If rho is not given as a mean field,
            # it will be automatically initialized as a function with all coefficients equal to zero. In that
            # case, a field with all coefficients equal to one is returned.
            # TODO: redo when restructuring the initialization process.
            if 'rho' in list(self._fieldDict.keys()) and sum(self._fieldDict['rho'].getCoefficientArray() ) != 0.:
                return self._fieldDict['rho'].getTensor()
            else:
                return self._oneField.getTensor()
        elif self.isMeanFlowVertexValuesClass():
            if 'rho' in list(self._fieldDict.keys()):
                return self._fieldDict['rho']
            else:
                return self._oneField

        else:
            if 'rho' in list(self._fieldDict.keys()):
                return self._fieldDict['rho']
            else:
                return self._zeroField

    @property
    def rhou(self):
        if 'rhou' in list(self._fieldDict.keys()):
            return self._fieldDict['rhou']
        else:
            return self._mean.rho * self.u + self.rho * self._mean.u

    def rhoY(self,species):
        return self._mean.rho * self.Y(species) + self.rho * self._mean.Y(species)

    @property
    def T(self):
        if self.isMeanFlowClass():
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T'].getTensor()
            else:
                return self._oneField
            
        elif self.isMeanFlowVertexValuesClass():
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T']
            else:
                return self._oneField

        else:
            if 'T' in list(self._fieldDict.keys()):
                return self._fieldDict['T']
            else:
                return self._zeroField

    @property
    def tau(self):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            mean_nu = self.nuTot
            mean_u = self.u
            tau_out = mean_nu * iGrad(mean_u)
            tau_out += iT(tau_out)
            # Sophie: this if-clause if not really necessary, in the incompressible case the term is just zero
            if not self._param.Case.SetOfEquations['Energy']['Equation'] == 'None':
                tau_out += -2.0/3.0 * mean_nu * \
                            iDiv(mean_u) * iIdentity(iGrad(mean_u))
                
        else:
            mean_nu = self._mean.nuTot
            mean_u = self._mean.u
            fluc_nu = self.nulam
            tau_out = mean_nu * iGrad(self.u) + \
                        fluc_nu * iGrad(mean_u)
            tau_out += iT(tau_out)
            # Sophie: this if-clause if not really necessary, in the incompressible case the term is just zero
            if not self._param.Case.SetOfEquations['Energy']['Equation'] == 'None':
                tau_out += -2.0/3.0 * mean_nu * iDiv(self.u) * iIdentity(iGrad(self.u))
                tau_out += -2.0/3.0 * fluc_nu * iDiv(mean_u) * iIdentity(iGrad(self.u))
        return tau_out

    @property
    def Tb(self):
        if self.isMeanFlowClass():
            if 'Tb' in list(self._fieldDict.keys()):
                return self._fieldDict['Tb'].getTensor()
        elif self.isMeanFlowVertexValuesClass():
            if 'Tb' in list(self._fieldDict.keys()):
                return self._fieldDict['Tb']
                
        else:
            return self._fieldDict['Tb']

    @property
    def Tu(self):
        if self.isMeanFlowClass():
            if 'Tu' in list(self._fieldDict.keys()):
                return self._fieldDict['Tu'].getTensor()
        elif self.isMeanFlowVertexValuesClass():
            if 'Tu' in list(self._fieldDict.keys()):
                return self._fieldDict['Tu']
                
        else:
            return self._fieldDict['Tu']

    @property
    def Tm(self):
        return self._fieldDict['Tm']

    @property
    def responseDomain(self):
        return self._fieldDict['responseDomain']

    @property
    def u(self):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            if 'u' in list(self._fieldDict.keys()):
                return self._fieldDict['u'].getTensor()
            else:
                return self._zeroVelocityField.getTensor()
        else:
            if 'u' in self._transportedQuantities:
                return  self._fieldDict['u']
            else:
                return self._zeroVelocityField

    @property
    def u_forcing_i(self):
        return self._fieldDict['u_forcing_i']

    @property
    def u_forcing_r(self):
        return self._fieldDict['u_forcing_r']

    @property
    def u_forcing(self):
        forcing = self.u_forcing_r + self.u_forcing_i * 1j
        return forcing.getTensor()

    @property
    def ut(self):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            if 'ut' in list(self._fieldDict.keys()):
                return self._fieldDict['u'][2]
            else:
                mesh = self._fieldDict[
                    list(self._fieldDict.keys())[0]].space.mesh
                return Constant(mesh, 0.0)
        else:
            if  self._param.Case.TransVelFluc:
                return self._fieldDict['u'][2]
            else:
                return Constant(0)

    def Y(self, specie):
        if self.isMeanFlowClass() or self.isMeanFlowVertexValuesClass():
            return self._fieldDict[specie].getTensor()
        else:
            return self._fieldDict[specie]
                        
