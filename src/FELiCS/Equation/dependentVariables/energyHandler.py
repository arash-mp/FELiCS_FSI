import numpy as np
from FELiCS.Equation.dependentVariables.viscosityModels import *

class energyHandler:
    """
    This class is used to link the temperature with the transported variable
    (e.g. enthalpy, sensible energy...). Depending on which state variable
    is linearized (or already known), the respective other is calculated.

    Parameters:
    ----------
    None

    Attributes:
    ----------
    None

    Methods:
    -------
    __init__():
        Initializing the class and link temperature with transported energy variable
        (e.g. enthalpy, sensible energy...)

    _relateConservativeToPrimitiveVariablesEnergy(mean='None'):
        Relates conservative variables to primitive variables for energy calculation.

    _getNeededFieldsForLinearEnergy():
        Returns the needed fields for linear energy calculation.

    _initializeMolecularHeatDiffusionFluctuation(mean='None'):
        Initializes molecular heat diffusion fluctuation.

    _additionalFieldsToBeReadEnergy():
        Returns additional fields to be read for energy calculation.

    """


    def __init__(self):
        """
        Initializing the class and link temperature with transported energy variable
        (e.g. enthalpy, sensible energy...)

        Parameters:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Returns:
        """
        pass

    def _relateConservativeToPrimitiveVariablesEnergy(
                                            self,
                                            mean='None',
                                            ):
        """
        Relates the conservative variables to primitive variables for energy equation.

        Parameters:
        - mean: FELiCS mean flow object

        Returns:
        """

        alreadyInitializedFields = list(self._fieldDict.keys())
        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']

        if mean == 'None':
            mean = self._mean

        ############################################### 
        ################### Enthalpy ##################
        ############################################### 
        if energyEquationType == 'Enthalpy':
            if sum(el in ['h', 'T'] for el in list(alreadyInitializedFields)) == 1:            
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
            if sum(el in ['progress', 'T'] for el in list(alreadyInitializedFields)) == 1:            
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
        
        ############################################### 
        ############### Energy-p linear ###############
        ############################################### 
        elif energyEquationType == 'primitive-p':
            # In this type of energy eq. we use only primitive 
            # variables so we don't need to define anything here.
            pass
        
        else: 
            raise Exception("Energy equation type " + energyEquationType + " not known.")

    def _getNeededFieldsForLinearEnergy(self):
        """
        Returns the list of fields needed for linear energy equation.

        Returns:
        - List of fields needed for linear energy equation.
        """

        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']
        if energyEquationType == 'Enthalpy':
            return ['T', 'h', 'alpha']
        if energyEquationType == 'primitive-p':
            return ['rho', 'p', 'T']
        elif energyEquationType == 'ProgressVariableLinear':
            return ['T', 'progress']

    def _initializeMolecularHeatDiffusionFluctuation(
                                    self,
                                    mean='None'
                                    ):
        """
        Initializes the molecular heat diffusion fluctuation.

        Parameters:
        - mean: FELiCS mean flow object

        Returns:
        """

        alreadyDeterminedFields = list(self._fieldDict.keys())
        viscosityModel = self._param.Case.MolViscPerturbModel
        if viscosityModel['type'] == 'Constant':
            pass
        elif viscosityModel['type'] == 'Sutherland mean':
            if 'rho' in alreadyDeterminedFields and not 'alpha' in alreadyDeterminedFields:
                if mean == 'None':
                    mean = self._mean
                if self._isSolution:
                    mean_alpha = mean.fieldDict['alpha']
                else:
                    mean_alpha = mean.alpha
                Ts = viscosityModel['Constants']['Ts']
                foobar, fluc = SutherlandFluctuationMean(mean, self.rho, Ts)
                self._fieldDict['alpha'] = mean_alpha * fluc

    def _additionalFieldsToBeReadEnergy(self):
        """
        Returns the list of additional fields to be read for energy equation.

        Returns:
        - List of additional fields to be read for energy equation.
        """

        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']
        if energyEquationType == 'Enthalpy': 
            return ['cp', 'alpha', 'he', 'T', 'molarMass']
        if energyEquationType == 'ProgressVariableLinear': 
            return ['T', 'Tu', 'Tb', 'rho']
        if energyEquationType == 'primitive-p': 
            # return ['rho', 'cp', 'T', 'p', 'gamma', 'Pr']
            return ['rho', 'cp', 'T', 'p', 'gamma']
        else:
            return []
