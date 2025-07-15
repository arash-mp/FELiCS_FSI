import FELiCS.Equation.Reactions.janafopenfoam as janafopenfoam


class Janafspecie:
    """
    Represents a chemical species using JANAF thermodynamic data.

    This class models the thermodynamic properties of a species based on 
    JANAF polynomial coefficients for specific heat capacity and entropy 
    over a range of temperatures.

    **Initialize the Janafspecie object**

    Parameters
    ----------
    W : float
        Molecular weight of the species [kg/kmol].
    Tlow : float
        Lower temperature limit for the thermodynamic data [K].
    Thigh : float
        Upper temperature limit for the thermodynamic data [K].
    Tcommon : float
        Common temperature separating the low and high temperature ranges [K].
    lowCpCoeffs : list of float
        Coefficients for the specific heat capacity polynomial in the low temperature range.
    highCpCoeffs : list of float
        Coefficients for the specific heat capacity polynomial in the high temperature range.
    P : float
        Reference pressure [Pa].

    Attributes
    ----------
    W : float
        Molecular weight [kg/kmol].
    R : float
        Specific gas constant [J/kg K].
    Tlow : float
        Lower temperature limit [K].
    Thigh : float
        Upper temperature limit [K].
    Tcommon : float
        Temperature boundary between low and high polynomial fits [K].
    lowCpCoeffs : list of float
        Low-temperature range Cp polynomial coefficients.
    highCpCoeffs : list of float
        High-temperature range Cp polynomial coefficients.
    s0 : float
        Standard entropy at reference state [J/kmol K].
    S0f : float
        Standard entropy at reference state per unit mass [J/kg K].
    """
    
    def __init__(self, W, Tlow, Thigh, Tcommon, lowCpCoeffs, highCpCoeffs, P):
        """
        Initializes the Janafspecie instance.

        Parameters
        ----------
        W : float
            Molecular weight of the species [kg/kmol].
        Tlow : float
            Lower temperature limit for the thermodynamic data [K].
        Thigh : float
            Upper temperature limit for the thermodynamic data [K].
        Tcommon : float
            Common temperature separating the low and high temperature ranges [K].
        lowCpCoeffs : list of float
            Coefficients for the specific heat capacity polynomial in the low temperature range.
        highCpCoeffs : list of float
            Coefficients for the specific heat capacity polynomial in the high temperature range.
        P : float
            Reference pressure [Pa].
        """


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
