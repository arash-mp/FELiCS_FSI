import FELiCS.Reactions.janafopenfoam as janafopenfoam


class Janafspecie:

    def __init__(self, W, Tlow, Thigh, Tcommon, lowCpCoeffs, highCpCoeffs, P):

        janaf = janafopenfoam.Janafopenfoam(P)

        self.W = W
        self.R = janaf.R_univ / self.W
        self.Tlow = Tlow
        self.Thigh = Thigh
        self.Tcommon = Tcommon
        self.highCpCoeffs = highCpCoeffs
        self.lowCpCoeffs = lowCpCoeffs

        self.s0 = janaf.janaf_s0(self)  # J/kmol K
        self.S0f = janaf.janaf_S0(self)  # J/kg K
