# Third party libraries
import numpy as np
from functions import printError, printDebug

class equationOfStateHandler:
    """
    This class is used to link the thermodynamic flow 
    variables (temperature, pressure, and density) 
    according to the ideal gas Equation of State (EoS).
    The class contains two versions of the EoS, a linearized
    one for the fluctuations, and a version for the mean
    flow variables.
    The user should give at least two of the three 
    thermodynamic variables: two mean fields in the mean
    flow file, and two transported variables for the 
    fluctuations.
    
    Currently, the following versions of the EoS are 
    implemented:
        .IdealGas
        .Low-Mach

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
        Initializing the class and link thermodynamic variables
        (pressure, density, temperature) for a perfect gas.

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        """
        pass
    
    def _initializeEoSFluctuations(self):
        '''
        This function is used to set the links 
        between the different fluctuations
        variables of the linearized equation of 
        state.
         
        Current types of linearized EoS:
            1. Low-Mach
            2. IdealGas (compressible)

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        '''
        # Fluctuations already initialized
        alreadyDefinedQuantities = list(self._fieldDict.keys())
        printDebug(True, f'-- Linearized EoS with input variables: {alreadyDefinedQuantities}')
        
        # Type of equation of state
        EoSType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if not EoSType in ['Low-Mach', 'IdealGas']:
            printError('Equation of State '+self._EoSModel()+\
                    ' not defined.')
            
        # Mean flow quantities available
        mean = self._mean

        # OPT.1 -- Calculate density from pressure and temperature
        if all(item in alreadyDefinedQuantities for item in ['p','T']):
            if self._isSolution:
                mean_T  = mean.T.x.array[:]
                mean_rho = mean.rho.x.array[:]
                if EoSType == 'IdealGas':
                    mean_Rspe = mean.R_spe.x.array[:]
            else:
                mean_T = mean.T
                mean_rho = mean.rho
                if EoSType == 'IdealGas':
                    mean_Rspe = mean.R_spe
                    
            if EoSType == 'Low-Mach':
                self._fieldDict['rho'] = -1 * mean_rho / mean_T * self.T
            elif EoSType == 'IdealGas':
                self._fieldDict['rho'] = (self.p - mean_rho*mean_Rspe*self.T)/(mean_Rspe*mean_T)

        
        # OPT.2 -- Calculate pressure from density and temperature
        elif all(item in alreadyDefinedQuantities for item in ['rho','T']):
            if EoSType == 'Low-Mach':
                printError('Equation of State '+self._EoSModel()+\
                    ' not defined to obtain p-fluctuations.')
            elif EoSType == 'IdealGas':
                if self._isSolution:
                    mean_T = mean.T.x.array[:]
                    mean_p = mean.p.x.array[:]
                    mean_Rspe = mean.R_spe.x.array[:]
                else:
                    mean_T = mean.T
                    mean_p = mean.p
                    mean_Rspe = mean.R_spe
                    
                self._fieldDict['p'] = mean_Rspe*(self.rho*mean_T + mean_p*self.T)
                    
        
        # OPT.3 -- Calculate temperature from density and pressure
        elif all(item in alreadyDefinedQuantities for item in ['rho','p']):
            if self._isSolution:
                mean_T = mean.fieldDict['T'].x.array[:]
                mean_rho = mean.fieldDict['rho'].x.array[:]
                # mean_T = mean.T.x.array[:]
                # mean_rho = mean.rho.x.array[:]
                if EoSType == 'IdealGas':
                    mean_Rspe = mean.fieldDict['R_spe'].x.array[:]
                    # mean_Rspe = mean.R_spe.x.array[:]
            else:
                mean_T = mean.T
                mean_rho = mean.rho
                if EoSType == 'IdealGas':
                    mean_Rspe = mean.R_spe
                    
            if EoSType == 'Low-Mach':
                self._fieldDict['T'] = -1 * self.rho / mean_rho * mean_T
            elif EoSType == 'IdealGas':
                self._fieldDict['T'] = (self.p - self.rho*mean_Rspe*mean_T)/(mean_Rspe*mean_rho)
                
    def _additionalFieldsToBeReadEoS(self):
        '''
        Define which additional fields must be read 
        from the mean flow file in order to apply the
        EoS fluctuation equation.
        '''
        EoSEquationType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if EoSEquationType == 'IdealGas': 
            return ['R_spe']
        else:
            return []

