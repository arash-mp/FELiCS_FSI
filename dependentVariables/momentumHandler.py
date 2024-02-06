# Third party libraries
import numpy as np

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

    def _getNeededFieldsLinear(self):
        return ['u', 'rhou','p']

    def _relateConservativeToPrimitiveVariables(self, mean = 'None'):
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

    def _additionalFieldsToBeReadEnergy(self):
        pass
