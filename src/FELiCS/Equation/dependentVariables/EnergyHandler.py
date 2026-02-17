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
# Third party libraries
import  numpy   as np

# Local Libraries and methods
from    FELiCS.Equation.dependentVariables.viscosityModels  import *

class EnergyHandler:
    """
    Handles the relationship between temperature and transported energy variables
    (e.g., enthalpy, sensible energy) for FELiCS simulations.

    Depending on which state variable is linearized or known, the respective other is calculated.

    **Initialize the energyHandler object**

    Parameters
    ----------
    None

    Attributes
    ----------
    _param : object
        FELiCS parameter object.
    _mean : object
        FELiCS mean flow object.
    _fieldDict : dict
        Dictionary of field variables.
    _isSolution : bool
        Flag indicating if the current field is a solution.
    """


    def __init__(self):
        """
        Initializes the energyHandler instance and links temperature with the transported energy variable.

        Parameters
        ----------
        None
        """
        pass

    def _relate_conservative_to_primitive_variables_energy(
        self,
        mean='None',
    ):
        """
        Relates conservative variables to primitive variables for energy calculation.

        Parameters
        ----------
        mean : object, optional
            FELiCS mean flow object. If not provided, uses self._mean.

        Raises
        ------
        Exception
            If required fields for enthalpy or progress variable calculation are missing.
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
                        mean_cp = mean.field_dict['cp']
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
                            mean_Tu  = mean.tu
                            mean_Tb  = mean.tb
                        else:
                            mean_Tu  = mean.tu
                            mean_Tb  = mean.tb
                        self._fieldDict['T'] = self.y('progress') * (mean_Tb - mean_Tu)
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

    def _get_needed_fields_for_linear_energy(self):
        """
        Returns the list of fields needed for linear energy equation.

        Returns
        -------
        fields : list of str
            List of fields needed for linear energy equation.
        """

        energyEquationType = self._param.Case.SetOfEquations['Energy']['Equation']
        if energyEquationType == 'Enthalpy':
            return ['T', 'h', 'alpha']
        if energyEquationType == 'primitive-p':
            return ['rho', 'p', 'T']
        elif energyEquationType == 'ProgressVariableLinear':
            return ['T', 'progress']

    def _initialize_molecular_heat_diffusion_fluctuation(
        self,
        mean='None',
    ):
        """
        Initializes the molecular heat diffusion fluctuation.

        Parameters
        ----------
        mean : object, optional
            FELiCS mean flow object. If not provided, uses self._mean.
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
                    mean_alpha = mean.field_dict['alpha']
                else:
                    mean_alpha = mean.alpha
                Ts = viscosityModel['Constants']['Ts']
                foobar, fluc = sutherland_fluctuation_mean(
                mean,
                self.rho,
                Ts,
                )
                self._fieldDict['alpha'] = mean_alpha * fluc

    def _additional_fields_to_be_read_energy(self):
        """
        Returns the list of additional fields to be read for energy equation.

        Returns
        -------
        fields : list of str
            List of additional fields to be read for energy equation.
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
