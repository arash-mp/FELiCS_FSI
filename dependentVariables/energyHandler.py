# Third party libraries
import numpy as np

class energyHandler:
    """
    This class is used to link the temperature with the transported variable
    (e.g. enthalpy, sensible energy...). Depending on which state variable
    is linearized (or alreaddy known), the respective other is calculated.

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
        """
        Initializing the class and link temperature with transported energy variable
        (e.g. enthalpy, sensible energy...)

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        """
        pass

    def _relateConservativeToPrimitiveVariablesEnergy(
                                            self,
                                            mean = 'None',
                                            ):
        alreadyInitializedFields = list(self._fieldDict.keys())
        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']

        if mean == 'None':
            mean = self._mean

        ############################################### 
        ################### Enthalpy ##################
        ############################################### 
        if energyEquationType == 'Enthalpy':
            if sum(el in ['h','T'] for el in list(alreadyInitializedFields)) == 1:            
                # so far cp is constant, which in reacting flows might be a strong assumption.
                if 'T' in alreadyInitializedFields:
                    if self._isSolution:
                        mean_cp = mean.fieldDict['cp']
                    else:
                        mean_cp = mean.cp
                    self._fieldDict['h'] = mean_cp * self._fieldDict['T']
                elif 'h' in alreadyInitializedFields:
                    raise Exception('Calculation of enthalpy, h, from temperature, T, not yet implemented. Check energyHandler.')
                else:
                    raise Exception('Enthalpy equation is chosen, however neither temperature, T, nor enthalpy, h, are available to calculate the respective other')

        ############################################### 
        ########### Progress Variable linear ##########
        ############################################### 
        elif energyEquationType == 'ProgressVariableLinear':
            if sum(el in ['progress','T'] for el in list(alreadyInitializedFields)) == 1:            
                    if 'progress' in alreadyInitializedFields:
                        if self._isSolution:
                            mean_Tu  = mean.Tu
                            mean_Tb  = mean.Tb
                        else:
                            mean_Tu  = mean.Tu
                            mean_Tb  = mean.Tb
                        self._fieldDict['T'] = self.Y('progress') * (mean_Tb - mean_Tu)
                    else:
                         raise Exception('ProgressVariableLinear is chosen, however, the progress variable is not available to calculate the temperature.')
        else: 
            raise Exception("Energy equation type " + energyEquationType + " not known.")

    def _getNeededFieldsForLinearEnergy(self):
        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']
        if energyEquationType == 'Enthalpy':
            return ['T', 'h']
        elif energyEquationType == 'ProgressVariableLinear':
            return ['T', 'progress']

    def _additionalFieldsToBeReadEnergy(self):
        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']
        if energyEquationType == 'Enthalpy': 
            return ['cp', 'alpha','he','T','molarMass']
        if energyEquationType == 'ProgressVariableLinear': 
            return ['T', 'Tu', 'Tb'] 
        if energyEquationType == 'primitive-p': 
            return ['cp', 'T', 'p', 'gamma', 'Pr'] 
        else:
            return []
