# Third party libraries
import numpy as np
from functions import printError

class equationOfStateHandler:
    """
    This class is used to provide the link between pressure, density and temperature fluctuation.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution

    Private attributes:

    Protected attributes:

    Public attributes:

    """

    def __init__(
            self,
            param,
            mean
    ):
        '''
        Adding pressure fluctuation to the fieldDict

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        '''
        
        alreadyDefinedQuantities = list(self._fieldDict.keys())

        print(alreadyDefinedQuantities)
        if all(item in alreadyDefinedQuantities for item in ['p','T']):
        # Calculate density from pressure and temperature
            if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != self.rho.shape[0]:
                local_mean_T = mean.T.compute_vertex_values()
                local_mean_rho = mean.rho.compute_vertex_values()
            else:
                local_mean_T = mean.T
                local_mean_rho = mean.rho
            if param.Case.SetOfEquations['EquationOfState']['Equation'] == 'Low-Mach':
                self._fieldDict['rho'] = -mean['rho']/mean['T']*self.T
            else:
                printError('Equation of State ' + param.Case.SetOfEquations['EquationOfState']['Equation'] + ' not defined.')

        elif all(item in alreadyDefinedQuantities for item in ['rho','T']):
        # Calculate pressure from density and temperature
            printError('Equation of state not defined to yield pressure fluctuation.')
            
        elif all(item in alreadyDefinedQuantities for item in ['rho','p']):
        # Calculate temperature from density and pressurea
            if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != self.rho.shape[0]:
                local_mean_T = mean.T.compute_vertex_values()
                local_mean_rho = mean.rho.compute_vertex_values()
            else:
                local_mean_T = mean.T
                local_mean_rho = mean.rho
            if param.Case.SetOfEquations['EquationOfState']['Equation'] == 'Low-Mach':
                self._fieldDict['T'] = -self.rho / local_mean_rho * local_mean_T
            else:
                printError('Equation of State ' + param.Case.SetOfEquations['EquationOfState']['Equation'] + ' not defined.')

