import numpy as np
from FELiCS.Equation.dependentVariables.viscosityModels import *

class momentumHandler:
    """
    Handles momentum-related calculations and variable transformations.

    This class provides methods for determining required fields, relating conservative
    and primitive variables, and initializing molecular momentum diffusion fluctuations
    based on the selected viscosity model.

    **Initialize the momentumHandler object**

    Parameters
    ----------
    None
    
    Attributes
    ----------
    _param : object
        Simulation parameter object containing case configuration.
    _fieldDict : dict
        Dictionary of available field variables.
    _mean : object
        Mean flow field object.
    _isSolution : bool
        Indicates if the handler is used for a solution field.
    """

    def __init__(self):
        """
        Initialize the momentumHandler instance.
        """
        pass

    def _getNeededFieldsForLinearMomentum(self):
        """
        Determine the fields required for linear momentum calculations.

        Returns
        -------
        list of str
            List of field names required for linear momentum calculations.
        """
        outList = ['u', 'rhou', 'p']
        if not self._param.Case.MolViscPerturbModel['type'] == 'Constant': 
            outList.append('nulam')
        return outList

    def _relateConservativeToPrimitiveVariablesMomentum(self, mean='None'):
        """
        Relate conservative variables to primitive variables for momentum calculations.

        Parameters
        ----------
        mean : str or object, optional
            Mean values to use in the calculations. If 'None', uses the default mean values.
        """
        alreadyDeterminedFields = list(self._fieldDict.keys())        
        if mean == 'None':
            mean = self._mean
        if sum(el in ['rhou', 'u'] for el in list(alreadyDeterminedFields)) == 1:
            if self._isSolution:
                mean_u = mean.fieldDict['u']
                if 'rho' in list(mean.fieldDict.keys()):
                    mean_rho = mean.fieldDict['rho']
                else:
                    mean_rho = mean.rho
            else:
                mean_u = mean.u
                mean_rho = mean.rho
            if 'rhou' in alreadyDeterminedFields:
                self._fieldDict['u'] = (self.rhou - mean_u * self.rho) / mean_rho
            elif 'u' in alreadyDeterminedFields:
                self._fieldDict['rhou'] = self.u * mean_rho + mean_u * self.rho

    def _initializeMolecularMomentumDiffusionFluctuation(self, mean='None'):
        """
        Initialize molecular momentum diffusion fluctuations based on the viscosity model.

        Parameters
        ----------
        mean : str or object, optional
            Mean values to use in the calculations. If 'None', uses the default mean values.

        Raises
        ------
        Exception
            If the viscosity model type is not implemented.
        """
        alreadyDeterminedFields = list(self._fieldDict.keys())
        viscosityModel = self._param.Case.MolViscPerturbModel
        if viscosityModel['type'] == 'Constant':
            pass
        elif viscosityModel['type'] == 'Sutherland mean':
            if 'rho' in alreadyDeterminedFields and not 'nulam' in alreadyDeterminedFields:
                if mean == 'None':
                    mean = self._mean
                mixture = self._param.Case.MolViscPerturbModel
                Ts = viscosityModel['Constants']['Ts']
                nulam, fluc = SutherlandFluctuationMean(mean, self.rho, Ts)
                self._fieldDict['nulam'] = nulam 
        else:
            raise Exception("Viscosity model " + viscosityModel['type'] + " not implemented.")

    def _additionalFieldsToBeReadEnergy(self):
        """
        Determine additional fields required for energy calculations.

        Notes
        -----
        This is a placeholder method and should be implemented as needed.
        """
        pass
