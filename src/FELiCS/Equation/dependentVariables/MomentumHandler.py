#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
import numpy as np
from FELiCS.Equation.dependentVariables.viscosityModels import *

class MomentumHandler:
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

    def _get_needed_fields_for_linear_momentum(self):
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

    def _relate_conservative_to_primitive_variables_momentum(self,
    mean='None',
    ):
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
                mean_u = mean.field_dict['u']
                if 'rho' in list(mean.field_dict.keys()):
                    mean_rho = mean.field_dict['rho']
                else:
                    mean_rho = mean.rho
            else:
                mean_u = mean.u
                mean_rho = mean.rho
            if 'rhou' in alreadyDeterminedFields:
                self._fieldDict['u'] = (self.rhou - mean_u * self.rho) / mean_rho
            elif 'u' in alreadyDeterminedFields:
                self._fieldDict['rhou'] = self.u * mean_rho + mean_u * self.rho

    def _initialize_molecular_momentum_diffusion_fluctuation(self,
    mean='None',
    ):
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
                nulam, fluc = sutherland_fluctuation_mean(
                mean,
                self.rho,
                Ts,
                )
                self._fieldDict['nulam'] = nulam 
        else:
            raise Exception("Viscosity model " + viscosityModel['type'] + " not implemented.")

    def _additional_fields_to_be_read_energy(self):
        """
        Determine additional fields required for energy calculations.

        Notes
        -----
        This is a placeholder method and should be implemented as needed.
        """
        pass
