class heatReleaseHandler:
    """
    This class is used to calculate the heat release.

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
            reaction
    ):
        """
        Adding heat release rate to the fieldDict

        Function arguments:
        - param: FELiCS reaction object

        Function returns:
        """
        Qtot, QList = reaction.dQ(self)
        self._fieldDict['Q'] = Qtot
        for Q in QList:
            self._fieldDict['Q' + str(QList.index(Q))] = Q