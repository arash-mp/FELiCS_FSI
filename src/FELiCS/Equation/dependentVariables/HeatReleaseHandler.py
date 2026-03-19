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
class HeatReleaseHandler:
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

    def __init__(
        self,
        reaction
    ):
        """
        Initialize the heatReleaseHandler instance and compute heat release rates.

        Parameters
        ----------
        reaction : object
            FELiCS reaction object used to compute heat release rates.
        """
        Qtot, QList = reaction.d_q(self)
        self._fieldDict['Q'] = Qtot
        for q in QList:
            self._fieldDict['Q' + str(QList.index(q))] = q
