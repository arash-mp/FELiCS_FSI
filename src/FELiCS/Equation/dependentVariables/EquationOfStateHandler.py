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
# Local Libraries and methods
from    FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class EquationOfStateHandler:
    """
    Handles the thermodynamic relationships between temperature, pressure,
    and density using different forms of the Equation of State (EoS),
    including linearized and mean-flow versions.

    This class links flow variables based on the chosen EoS, which may
    be either `IdealGas` or `Low-Mach`. It ensures that, given any two
    of the three thermodynamic quantities (pressure, temperature, and
    density), the third can be deduced according to the selected EoS.

    It is primarily intended to support FELiCS simulations by connecting
    transported fluctuation variables with background mean fields.

    **Initialize the equationOfStateHandler object**

    Parameters
    ----------
    None

    Attributes
    ----------
    _fieldDict : dict
        Stores fluctuation variables computed from EoS.
    _param : object
        Reference to the FELiCS parameter configuration.
    _mean : object
        Reference to the FELiCS mean flow object.
    _isSolution : bool
        Indicates if the current instance deals with a fluctuation solution.
    """


    def __init__(self):
        """
        Initializes the equationOfStateHandler instance.

        Sets up references to FELiCS parameter and mean flow objects,
        and prepares internal data structures for managing fluctuation
        and mean field variables.
        """

        pass
    
    def _initialize_eo_s_fluctuations(
        self,
        mean = 'None',
    ):
        """
        Initializes the linearized Equation of State for fluctuations.

        Based on the given mean fields and already defined fluctuation
        variables, computes the missing thermodynamic variable using
        the specified EoS. Requires that exactly two of the three
        (p, T, rho) are available.

        Parameters
        ----------
        mean : object, optional
            FELiCS mean flow object. If not provided, the default internal
            `_mean` object is used.

        Raises
        ------
        Exception
            If an unsupported EoS is specified or insufficient variables
            are available to compute the third thermodynamic quantity.
        """

        # Fluctuations already initialized
        alreadyDefinedQuantities = list(self._fieldDict.keys())

        #only do something if two out of the three variables (p,T,rho) are aleardy determined
        if sum(el in ['rho','p', 'T'] for el in list(alreadyDefinedQuantities)) == 2:
            
            # Type of equation of state
            EoSType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
            if EoSType not in ['Low-Mach', 'IdealGas']:
                logger.error('Equation of State '+self._EoSType+' not defined.')
                raise Exception('Equation of State '+self._EoSType+' not defined.')
                
            # Get the mean flow class     
            if mean == 'None':
                mean = self._mean

            # OPT.1 -- Calculate density from pressure and temperature
            if all(item in alreadyDefinedQuantities for item in ['p','T']):
                logger.debug('Linearized EoS with input variables: [p, T] -> rho')
                if self._isSolution:
                    mean_T  = mean.field_dict['T']
                    mean_rho = mean.field_dict['rho']
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.field_dict['R_spe']
                else:
                    mean_T = mean.T
                    mean_rho = mean.rho
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.r_spe
                        
                if EoSType == 'Low-Mach':
                    self._fieldDict['rho'] = -1 * mean_rho / mean_T * self.T
                elif EoSType == 'IdealGas':
                    self._fieldDict['rho'] = (self.p - mean_rho*mean_Rspe*self.T)/(mean_Rspe*mean_T)

            
            # OPT.2 -- Calculate pressure from density and temperature
            elif all(item in alreadyDefinedQuantities for item in ['rho','T']):
                logger.debug('Linearized EoS with input variables: [rho, T] -> p')
                if EoSType == 'Low-Mach':
                    logger.error('Equation of State '+self._EoSModel()+' not defined to obtain p-fluctuations.')
                    raise Exception('Equation of State '+self._EoSModel()+' not defined to obtain p-fluctuations.')
                elif EoSType == 'IdealGas':
                    if self._isSolution:
                        mean_T = mean.field_dict['T']
                        mean_rho = mean.field_dict['rho']
                        mean_Rspe = mean.field_dict['R_spe']
                    else:
                        mean_T = mean.T
                        mean_rho = mean.rho
                        mean_Rspe = mean.r_spe
                        
                    self._fieldDict['p'] = mean_Rspe*(self.rho*mean_T + mean_rho*self.T)
                        
            
            # OPT.3 -- Calculate temperature from density and pressure
            elif all(item in alreadyDefinedQuantities for item in ['rho','p']):
                logger.debug('Linearized EoS with input variables: [rho, p] -> T')
                if self._isSolution:
                    mean_T = mean.field_dict['T']
                    mean_rho = mean.field_dict['rho']
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.field_dict['R_spe']
                else:
                    mean_T = mean.T
                    mean_rho = mean.rho
                    if EoSType == 'IdealGas':
                        mean_Rspe = mean.r_spe
                        
                if EoSType == 'Low-Mach':
                    self._fieldDict['T'] = -1 * self.rho / mean_rho * mean_T
                elif EoSType == 'IdealGas':
                    self._fieldDict['T'] = (self.p - self.rho*mean_Rspe*mean_T)/(mean_Rspe*mean_rho)
                    
                    
    def _get_needed_fields_for_linear_eo_s(self):
        """
        Returns the list of thermodynamic fields needed for linearized EoS.

        Based on the EoS type, returns the minimal set of variables required
        to apply the linearized equation to fluctuations.

        Returns
        -------
        list of str
            List of field names such as ['rho', 'p', 'T'].

        Raises
        ------
        Exception
            If the EoS type is not recognized.
        """

        EoSEquationType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if EoSEquationType == 'Low-Mach':
            return ['rho','T']
        if EoSEquationType == 'IdealGas':
            return ['rho', 'p', 'T']
        if EoSEquationType == 'None':
            return []
        else:
            raise Exception('Equation of state type ' + EoSEquationType + ' not implemented.')

                
    def _additional_fields_to_be_read_eo_s(self):
        """
        Specifies the additional mean flow fields required for EoS evaluation.

        Depending on the EoS model, returns a list of extra variables needed
        from the mean flow file (e.g., specific gas constant).

        Returns
        -------
        list of str
            List of additional required mean flow fields.

        Raises
        ------
        Exception
            If the EoS type is not implemented.
        """
        EoSEquationType = self._param.Case.SetOfEquations['EquationOfState']['Equation']
        if EoSEquationType == 'IdealGas': 
            return ['R_spe']
        elif EoSEquationType == 'Low-Mach':
            return ['rho','T']
        elif EoSEquationType == 'None':
            return []
        else:
            raise Exception('Equation of state type ' + EoSEquationType + ' not implemented.')

