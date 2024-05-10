# Third party libraries
import numpy as np
from FELiCS.Equation.dependentVariables.viscosityModels import *

class momentumHandler:
    """
    This class is used to manage the variables for the linear momentum 
    equation. Depending on which state variables are linearized (or 
    alreaddy known), the respective others are calculated.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution

    Private attributes:

    Protected attributes:

    Public attributes:

    """


    def __init__(self):
        """
        Initializing the class and link temperature with transported energy variable
        (e.g. enthalpy, sensible energy...)

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        """
        pass

    def _getNeededFieldsForLinearMomentum(self):
        outList = ['u', 'rhou','p']
        if not self._param.Case.Mixture.Viscosity['type'] == 'Constant': 
            outList.append('nulam')
        return outList

    def _relateConservativeToPrimitiveVariablesMomentum(self, mean = 'None'):
        '''
        This function checks if conservative and/or primitive variable is already defined. 
        If only one of them is, it relates the one to the respective other.
        '''
        alreadyDeterminedFields = list(self._fieldDict.keys())        
        if mean == 'None':
            mean = self._mean
        # check if EITHER rhou OR u is already defined if not calculate the other
        if sum(el in ['rhou','u'] for el in list(alreadyDeterminedFields)) == 1:
            #if self is a solution, then the coefficients are used for the calculation...
            if self._isSolution:
                mean_u = mean.fieldDict['u']
                # The following if clause is a quick fix. Should be removed once the mean flow class has been changed to Tensor notation
                if 'rho' in list(mean.fieldDict.keys()):
                    mean_rho = mean.fieldDict['rho']
                else:
                    #mean_rho = mean.oneField
                    mean_rho = mean.rho
            # else it is a tensor object and will be used to build a ufl formulation
            else:
                mean_u = mean.u
                mean_rho = mean.rho
            # Calculate the primitive OR the conservative variable from the respective other
            if 'rhou' in alreadyDeterminedFields:
                self._fieldDict['u'] = (self.rhou - mean_u * self.rho) /  mean_rho
            elif 'u' in alreadyDeterminedFields:
                self._fieldDict['rhou'] = self.u * mean_rho + mean_u * self.rho

    def _initializeMolecularMomentumDiffusionFluctuation(
                                    self,
                                    mean = 'None'
                                    ):
        alreadyDeterminedFields = list(self._fieldDict.keys())
        viscosityModel = self._param.Case.Mixture.Viscosity
        if viscosityModel['type'] == 'Constant':
            pass
        elif viscosityModel['type'] == 'Sutherland mean':
            if 'rho' in alreadyDeterminedFields and not 'nulam' in alreadyDeterminedFields:
                if mean == 'None':
                    mean = self._mean
                Mixture = self._param.Case.Mixture
                Ts = Mixture.Viscosity['Constants']['Ts']
                nulam,fluc = SutherlandFluctuationMean(mean, self.rho, Ts)
                self._fieldDict['nulam'] = nulam 
        else:
            raise Exception("Viscosity model " + viscosityModel['type'] + " not implemented.")

    def _additionalFieldsToBeReadEnergy(self):
        pass
