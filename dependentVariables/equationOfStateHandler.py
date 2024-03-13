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
    
    def _initializeEoSFluctuations(
                            self,
                            mean = 'None'
                            ):
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

        #only do something if two out of the three variables (p,T,rho) are aleardy determined
        if sum(el in ['rho','p', 'T'] for el in list(alreadyDefinedQuantities)) == 2:
            
            # Type of equation of state
            EoSType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
            if not EoSType in ['Low-Mach', 'IdealGas']:
                printError('Equation of State '+self._EoSType+\
                        ' not defined.')
                
            # Get the mean flow class     
            if mean == 'None':
                mean = self._mean

            # OPT.1 -- Calculate density from pressure and temperature
            if all(item in alreadyDefinedQuantities for item in ['p','T']):
                printDebug(True, '-- Linearized EoS with input variables: [p, T] -> rho')
                if self._isSolution:
                    mean_T  = mean.fieldDict['T']
                    mean_rho = mean.fieldDict['rho']
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.fieldDict['R_spe']
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
                printDebug(True, '-- Linearized EoS with input variables: [rho, T] -> p')
                if EoSType == 'Low-Mach':
                    printError('Equation of State '+self._EoSModel()+\
                        ' not defined to obtain p-fluctuations.')
                elif EoSType == 'IdealGas':
                    if self._isSolution:
                        mean_T = mean.fieldDict['T']
                        mean_p = mean.fieldDict['p']
                        mean_Rspe = mean.fieldDict['R_spe']
                    else:
                        mean_T = mean.T
                        mean_p = mean.p
                        mean_Rspe = mean.R_spe
                        
                    self._fieldDict['p'] = mean_Rspe*(self.rho*mean_T + mean_p*self.T)
                        
            
            # OPT.3 -- Calculate temperature from density and pressure
            elif all(item in alreadyDefinedQuantities for item in ['rho','p']):
                printDebug(True, '-- Linearized EoS with input variables: [rho, p] -> T')
                if self._isSolution:
                    mean_T = mean.fieldDict['T']
                    mean_rho = mean.fieldDict['rho']
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.fieldDict['R_spe']
                else:
                    mean_T = mean.T
                    mean_rho = mean.rho
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.R_spe
                        
                if EoSType == 'Low-Mach':
                    self._fieldDict['T'] = -1 * self.rho / mean_rho * mean_T
                elif EoSType == 'IdealGas':
                    self._fieldDict['T'] = (self.p - self.rho*mean_Rspe*mean_T)/(mean_Rspe*mean_rho)
                    
                    
    def _getNeededFieldsForLinearEoS(self):
        '''This is a standard function for handlers, which defines the additional Fields necessary to be determinied. 
        The EoS only clsoses the variables appearing in the other equations. So it is not necessary to determine
        additional fields. Therefore, this function only returns an empty list, but still exists and is called for 
        reasons of consistency'''
        EoSEquationType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if EoSEquationType == 'Low-Mach':
            return ['rho','T']
        if EoSEquationType == 'IdealGas':
            return ['rho', 'p', 'T']
        if EoSEquationType == 'None':
            return []
        else:
            raise Exception('Equation of state type ' + EoSEquationType + ' not implemented.')

                
    def _additionalFieldsToBeReadEoS(self):
        '''
        Define which additional fields must be read 
        from the mean flow file in order to apply the
        EoS fluctuation equation.
        '''
        EoSEquationType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if EoSEquationType == 'IdealGas': 
            return ['R_spe']
        elif EoSEquationType == 'Low-Mach':
            return ['rho','T']
        elif EoSEquationType == 'None':
            return []
        else:
            raise Exception('Equation of state type ' + EoSEquationType + ' not implemented.')

