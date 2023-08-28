class temperatureHandler:
    """
    This class is used to calculate the temperature.

    Parent classes:

    Child classes:
    - fluctuationClass
    - fluctuationSolution

    Private attributes:

    Protected attributes:

    Public attributes:

    """


    def __init__(self, param, mean):
        """
        Adding temperature to the fieldDict

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        """

        self._fieldDict['T'] = self.getTemp(
            param,
            mean,
            self._meanfieldDict if hasattr(self, "_meanfieldDict") else None,
            self.rho,
            self.p
        )

    def getTemp(self, param, mean, ClassDict, rho, p):
        import numpy as np

        if param.Case.Reaction and not \
        param.Case.Mixture.getReactionMechanism()['type'] in ['2S-SM2']:
            # These hardcoded values are only temporary/unfinished
            if hasattr(mean.molarMass, 'vector'):
                temp = mean.molarMass.vector[0]
            else:
                temp = mean.molarMass[0]
            if temp == 0:
                M = 28.949
                print('using constant molarMass')
            else:
                M = mean.molarMass * 1e3
                print('using calculated molarMass field')
            pRef = 101300
            R = 8314.4598 / M

            if isinstance(mean.rho, np.ndarray) and mean.rho.shape[0] != \
                    rho.shape[0]:
                local_rho = ClassDict['rho'].compute_vertex_values()
            else:
                local_rho = mean.rho
            Tm = pRef / R / local_rho
            flucT = -rho / local_rho * Tm

        else:
            if isinstance(mean.rho, np.ndarray) and mean.rho.shape[0] != \
                    rho.shape[0]:
                local_rho = ClassDict['rho'].compute_vertex_values()
                local_T = ClassDict['T'].compute_vertex_values()
                if param.Case.Compressible:
                    local_p = ClassDict['p'].compute_vertex_values()
            else:
                local_rho = mean.rho
                local_T = mean.T
                if param.Case.Compressible:
                    local_p = mean.p
            flucT = -rho / local_rho * local_T
            if param.Case.Compressible:
                flucT += p / local_p * local_T

        return flucT
