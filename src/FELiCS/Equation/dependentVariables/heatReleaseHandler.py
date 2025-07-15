class heatReleaseHandler:
    """
    Calculates the heat release for combustion simulations.

    This class computes the total and component-wise heat release rates
    and stores them in the internal field dictionary. It is typically used as a
    component of `fluctuationClass` and `fluctuationSolution`.

    **Initialize the heatReleaseHandler object**

    Parameters
    ----------
    reaction : object
        FELiCS reaction object used to compute heat release rates.

    Attributes
    ----------
    _fieldDict : dict
        Dictionary storing total and component heat release rates.
    """

    def __init__(self, reaction):
        """
        Initialize the heatReleaseHandler instance and compute heat release rates.

        Parameters
        ----------
        reaction : object
            FELiCS reaction object used to compute heat release rates.
        """
        Qtot, QList = reaction.dQ(self)
        self._fieldDict['Q'] = Qtot
        for Q in QList:
            self._fieldDict['Q' + str(QList.index(Q))] = Q