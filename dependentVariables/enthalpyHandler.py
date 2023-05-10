# Third party libraries
import numpy as np


class enthalpyHandler:
    """
    This class provides functions for calculating the enthalpy.

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
        """
        Adding enthalpy to the fieldDict

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object

        Function returns:
        """

        if not param.Case.Mixture.ReactionMechanism['type'] == 'NOx':
            if isinstance(mean.cp, np.ndarray):
                if mean.cp.shape[0] != self._fieldDict['T'].shape[0]:
                    self._fieldDict['h'] \
                        = self._meanfieldDict['cp'].compute_vertex_values() \
                          * self._fieldDict['T']
                else:
                    self._fieldDict['h'] = mean.cp * self._fieldDict['T']
            else:
                self._fieldDict['h'] = mean.cp * self._fieldDict['T']