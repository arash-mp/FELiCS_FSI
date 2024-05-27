import dolfin as do

# noinspection PyUnresolvedReferences
import FELiCS.Equation.Reactions.janafspecie as janafspecie
# noinspection PyUnresolvedReferences
import FELiCS.Equation.Reactions.janafopenfoam as janafopenfoam

"""
The mechanism is actually called 2S-CM2, but we should avoid python modules
that start with a number or contain dashes, so c2sm2 it is...
"""

class C2SM2:
    """
    C2SM2 class representing a chemical reaction mechanism.

    Attributes
    ----------
    YCH4_lim : float
        The threshold value for YCH4.
    P : int
        The order of basis functions.
    janaf : Janafopenfoam
        Instance of the Janafopenfoam class.
    N2, CO2, O2, CH4, CO, H2O, AR : Janafspecie
        Instances of the Janafspecie class representing different species.
    h0r_1, h0r_2 : float
        Enthalpy of reaction 1 and 2.
    s0r_1, s0r_2 : float
        Entropy of reaction 1 and 2.
    A1, A2, Ta1, Ta2 : float
        Reaction constants.
    pa : float
        Pressure in N/m^2.
    R, R_mol : float
        Universal gas constants.
    n_CH4_1, n_O2_1, n_CO_1, n_H2O_1 : float
        Stoichiometric coefficients for reaction 1.
    n_CO_2, n_O2_2, n_CO2_2 : float
        Stoichiometric coefficients for reaction 2.
    nu_CH4_1, nu_O2_1, nu_CO_1, nu_H2O_1 : float
        Kinematic exponents for reaction 1.
    nu_CO_2, nu_O2_2, nu_CO2_2 : float
        Kinematic exponents for reaction 2.
    Q1dT, Q1drho, Q1dYCH4, Q1dYO2, Q2dT, Q2drho, Q2dYCO, Q2dYO2, Q2dYCO2 : float
        Sensitivity variables.
    """

    def __init__(self, YCH4_lim, P):
        """
        Initialize the C2SM2 class.

        Parameters
        ----------
        YCH4_lim : float
            The threshold value for YCH4.
        P : int
            The order of basis functions.

        Returns
        -------
        None
        """
        self.janaf = janafopenfoam.Janafopenfoam(P)
        self.YCH4_lim = YCH4_lim
        self.P = P

        N2 = janafspecie.Janafspecie(28.0134, 200, 5000, 1000,
                                     [3.29868, 0.00140824, -3.96322e-06, 5.64152e-09, -2.44486e-12, -1020.9, 3.95037],
                                     [2.92664, 0.00148798, -5.68476e-07, 1.0097e-10, -6.75335e-15, -922.798, 5.98053], P)
        CO2 = janafspecie.Janafspecie(44.01, 200, 5000, 1000,
                                      [2.27572, 0.00992207, -1.04091e-05, 6.86669e-09, -2.11728e-12, -48373.1, 10.1885],
                                      [4.45362, 0.00314017, -1.27841e-06, 2.394e-10, -1.66903e-14, -48967, -0.955396], P)
        O2 = janafspecie.Janafspecie(31.9988, 200, 5000, 1000,
                                     [3.21294, 0.00112749, -5.75615e-07, 1.31388e-09, -8.76855e-13, -1005.25, 6.03474],
                                     [3.69758, 0.00061352, -1.25884e-07, 1.77528e-11, -1.13644e-15, -1233.93, 3.18917], P)
        CH4 = janafspecie.Janafspecie(16.043, 200, 5000, 1000,
                                      [0.778741, 0.0174767, -2.78341e-05, 3.04971e-08, -1.22393e-11, -9825.23, 13.7222],
                                      [1.68348, 0.0102372, -3.87513e-06, 6.78559e-10, -4.50342e-14, -10080.8, 9.6234], P)
        CO = janafspecie.Janafspecie(28.0106, 200, 5000, 1000,
                                     [3.26245, 0.00151194, -3.88176e-06, 5.58194e-09, -2.47495e-12, -14310.5, 4.8489],
                                     [3.02508, 0.00144269, -5.63083e-07, 1.01858e-10, -6.91095e-15, -14268.4, 6.10822], P)
        H2O = janafspecie.Janafspecie(18.0153, 200, 5000, 1000,
                                      [3.38684, 0.00347498, -6.3547e-06, 6.96858e-09, -2.50659e-12, -30208.1, 2.59023],
                                      [2.67215, 0.00305629, -8.73026e-07, 1.201e-10, -6.39162e-15, -29899.2, 6.86282], P)
        AR = janafspecie.Janafspecie(39.948, 200, 5000, 1000,
                                     [2.5, 0, 0, 0, 0, -745.375, 4.366],
                                     [2.5, 0, 0, 0, 0, -745.375, 4.366], P)

        self.N2 = N2
        self.CO2 = CO2
        self.O2 = O2
        self.CH4 = CH4
        self.CO = CO
        self.H2O = H2O
        self.AR = AR

        self.n_CH4_1 = 1
        self.n_O2_1 = 1.5
        self.n_CO_1 = 1
        self.n_H2O_1 = 2

        self.n_CO_2 = 1
        self.n_O2_2 = 0.5
        self.n_CO2_2 = 1

        self.nu_CH4_1 = 0.9
        self.nu_O2_1 = 1.1
        self.nu_CO_1 = 1
        self.nu_H2O_1 = 1

        self.nu_CO_2 = 1
        self.nu_O2_2 = 0.5
        self.nu_CO2_2 = 1

        self.A1 = 2e12
        self.Ta1 = 17611.7

        self.A2 = 6.32456e7
        self.Ta2 = 6038.29

        self.pa = 1e5
        self.R = self.janaf.R_univ
        self.R_mol = 8.3144598

        self.h0r_1 = -self.n_CH4_1 * self.janaf.janaf_hc(self.CH4) - self.n_O2_1 * self.janaf.janaf_hc(self.O2) + self.n_CO_1 * self.janaf.janaf_hc(self.CO) + self.n_H2O_1 * self.janaf.janaf_hc(self.H2O)
        self.h0r_2 = -self.n_CO_2 * self.janaf.janaf_hc(self.CO) - self.n_O2_2 * self.janaf.janaf_hc(self.O2) + self.n_CO2_2 * self.janaf.janaf_hc(self.CO2)

        self.s0r_1 = -self.n_CH4_1 * self.CH4.s0 - self.n_O2_1 * self.O2.s0 + self.n_CO_1 * self.CO.s0 + self.n_H2O_1 * self.H2O.s0
        self.s0r_2 = -self.n_CO_2 * self.CO.s0 - self.n_O2_2 * self.O2.s0 + self.n_CO2_2 * self.CO2.s0

        self.Q1dT = None
        self.Q1drho = None
        self.Q1dYCH4 = None
        self.Q1dYO2 = None
        self.Q2dT = None
        self.Q2drho = None
        self.Q2dYCO = None
        self.Q2dYO2 = None
        self.Q2dYCO2 = None

    def computeSensitivities(self, T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2):
        """
        Compute all the sensitivities of the progress rates and store them.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        bool
            Returns True when sensitivities are computed.
        """
        self.Q1dT = self.Q1dT_(T, rho, Y_CH4, Y_O2)
        self.Q1drho = self.Q1drho_(T, rho, Y_CH4, Y_O2)
        self.Q1dYCH4 = self.Q1dYCH4_(T, rho, Y_CH4, Y_O2)
        self.Q1dYO2 = self.Q1dYO2_(T, rho, Y_CH4, Y_O2)
        self.Q2dT = self.Q2dT_(T, rho, Y_CO, Y_O2, Y_CO2)
        self.Q2drho = self.Q2drho_(T, rho, Y_CO, Y_O2, Y_CO2)
        self.Q2dYCO = self.Q2dYCO_(T, rho, Y_O2)
        self.Q2dYO2 = self.Q2dYO2_(T, rho, Y_CO, Y_O2)
        self.Q2dYCO2 = self.Q2dYCO2_(T, rho)

        return True

    def deltaGT02_(self, T):
        """
        Compute the standard-state Gibbs function change for reversible reaction 2.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the Gibbs function change.
        """
        deltaGT02 = do.Expression("-n_O2_2*g_O2 - n_CO_2*g_CO + n_CO2_2*g_CO2",
                                  n_O2_2=self.n_O2_2, g_O2=self.janaf.janaf_g_expr(self.O2, T),
                                  n_CO_2=self.n_CO_2, g_CO=self.janaf.janaf_g_expr(self.CO, T),
                                  n_CO2_2=self.n_CO2_2, g_CO2=self.janaf.janaf_g_expr(self.CO2, T),
                                  degree=self.P)
        return deltaGT02

    def Kf1_(self, T):
        """
        Compute the forward rate of reaction 1.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the forward rate of reaction 1.
        """
        Kf1 = do.Expression("A1*exp(-Ta1/T)", A1=self.A1, Ta1=self.Ta1, T=T, degree=self.P)
        return Kf1

    def Kf2_(self, T):
        """
        Compute the forward rate of reaction 2.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the forward rate of reaction 2.
        """
        Kf2 = do.Expression("A2*exp(-Ta2/T)", A2=self.A2, Ta2=self.Ta2, T=T, degree=self.P)
        return Kf2

    def Kp2_(self, T):
        """
        Compute the equilibrium constant for reaction 2.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the equilibrium constant for reaction 2.
        """
        Kp2 = do.Expression("exp(-deltaGT0/(R*T))", deltaGT0=self.deltaGT02_(T), R=self.R, T=T, degree=self.P)
        return Kp2

    def Kc2_(self, T):
        """
        Compute the concentration equilibrium constant for reaction 2.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the concentration equilibrium constant for reaction 2.
        """
        Kc2 = do.Expression("pow((pa/(R*T)), (n_CO2_2-n_O2_2-n_CO_2)) * Kp2", pa=self.pa, R=self.R, T=T,
                            n_CO2_2=self.n_CO2_2, n_O2_2=self.n_O2_2, n_CO_2=self.n_CO_2,
                            Kp2=self.Kp2_(T), degree=self.P)
        return Kc2

    def Kr2_(self, T):
        """
        Compute the reverse rate of reaction 2.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the reverse rate of reaction 2.
        """
        Kr2 = do.Expression("Kf2 / Kc2", Kf2=self.Kf2_(T), Kc2=self.Kc2_(T), degree=self.P)
        return Kr2

    def Q1_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the rate of progress of reaction 1.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the rate of progress of reaction 1.
        """
        Q1 = do.Expression(
            "Kf1 * pow(rho_, nu_CH4_1 + nu_O2_1) * pow(Y_CH4_ / W_CH4, nu_CH4_1) * pow(Y_O2_ / W_O2, nu_O2_1)",
            Kf1=self.Kf1_(T), rho_=rho, Y_CH4_=Y_CH4, Y_O2_=Y_O2,
            nu_CH4_1=self.nu_CH4_1, nu_O2_1=self.nu_O2_1, W_CH4=self.CH4.W, W_O2=self.O2.W,
            degree=self.P)
        return Q1

    def Q2f_(self, T, rho, Y_CO, Y_O2):
        """
        Compute the forward rate of progress of reaction 2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the forward rate of progress of reaction 2.
        """
        Q2f = do.Expression(
            "Kf2 * pow(rho_, nu_CO_2 + nu_O2_2) * pow(Y_CO_ / W_CO, nu_CO_2) * pow(Y_O2_ / W_O2, nu_O2_2)",
            Kf2=self.Kf2_(T), rho_=rho, Y_CO_=Y_CO, Y_O2_=Y_O2,
            nu_CO_2=self.nu_CO_2, nu_O2_2=self.nu_O2_2, W_CO=self.CO.W, W_O2=self.O2.W,
            degree=self.P)
        return Q2f

    def Q2r_(self, T, rho, Y_CO2):
        """
        Compute the reverse rate of progress of reaction 2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the reverse rate of progress of reaction 2.
        """
        Q2r = do.Expression("Kr2 * pow(rho_, nu_CO2_2) * pow(Y_CO2_ / W_CO2, nu_CO2_2)",
                            Kr2=self.Kr2_(T), rho_=rho, Y_CO2_=Y_CO2,
                            nu_CO2_2=self.nu_CO2_2, W_CO2=self.CO2.W,
                            degree=self.P)
        return Q2r

    def Q2_(self, T, rho, Y_CO, Y_O2, Y_CO2):
        """
        Compute the rate of progress of reaction 2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the rate of progress of reaction 2.
        """
        Q2 = do.Expression("Q2f - Q2r",
                           Q2f=self.Q2f_(T, rho, Y_CO, Y_O2), Q2r=self.Q2r_(T, rho, Y_CO2), degree=self.P)
        return Q2

    def Q_(self, T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2):
        """
        Compute the heat-release rate of both reactions combined.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the combined heat-release rate.
        """
        dQ = do.Expression("-(Q1 * h0r_1 + Q2 * h0r_2)",
                           Q1=self.Q1_(T, rho, Y_CH4, Y_O2), Q2=self.Q2_(T, rho, Y_CO, Y_O2, Y_CO2),
                           h0r_1=self.h0r_1, h0r_2=self.h0r_2, degree=self.P)
        return dQ

    def Q(self, T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2):
        """
        Return the heat release rate as a Function.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Function
            Projected heat release rate function.
        """
        return do.project(self.Q_(T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2), Y_CH4.function_space())

    def dQ(self, fluc):
        """
        Return the fluctuations in heat release.

        Parameters
        ----------
        fluc : object
            Object containing fluctuation data for temperature, density, and species mass fractions.

        Returns
        -------
        tuple
            A tuple containing the total fluctuation in heat release and a list of individual fluctuations.
        """
        QDict = {}
        dQ1 = - self.h0r_1 * (self.Q1dT * fluc.T + self.Q1drho * fluc.rho + self.Q1dYCH4 * fluc.Y('CH4') + self.Q1dYO2 * fluc.Y('O2'))
        dQ2 = - self.h0r_2 * (self.Q2dT * fluc.T + self.Q2drho * fluc.rho + self.Q2dYCO * fluc.Y('CO') + self.Q2dYCO2 * fluc.Y('CO2') + self.Q2dYO2 * fluc.Y('O2'))
        dQ = dQ1 + dQ2
        return dQ, [dQ1, dQ2]

    def rrCH4_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the reaction rate of CH4.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the reaction rate of CH4.
        """
        rrCH4 = do.Expression("-n_CH4_1 * W_CH4 * Q1",
                              n_CH4_1=self.n_CH4_1, W_CH4=self.CH4.W,
                              Q1=self.Q1_(T, rho, Y_CH4, Y_O2), degree=self.P)
        return rrCH4

    def rrH2O_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the reaction rate of H2O.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the reaction rate of H2O.
        """
        rrH2O = do.Expression("n_H2O_1 * W_H2O * Q1",
                              n_H2O_1=self.n_H2O_1, W_H2O=self.H2O.W,
                              Q1=self.Q1_(T, rho, Y_CH4, Y_O2), degree=self.P)
        return rrH2O

    def rrO2_(self, T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2):
        """
        Compute the reaction rate of O2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the reaction rate of O2.
        """
        rrO2 = do.Expression("W_O2 * (- n_O2_1 * Q1 - n_O2_2 * Q2)",
                             n_O2_1=self.n_O2_1, n_O2_2=self.n_O2_2, W_O2=self.O2.W,
                             Q1=self.Q1_(T, rho, Y_CH4, Y_O2), Q2=self.Q2_(T, rho, Y_CO, Y_O2, Y_CO2), degree=self.P)
        return rrO2

    def rrCO_(self, T, rho, Y_CH4, Y_CO, Y_O2, Y_CO2):
        """
        Compute the reaction rate of CO.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the reaction rate of CO.
        """
        rrCO = do.Expression("W_CO * (n_CO_1 * Q1 - n_CO_2 * Q2)",
                             n_CO_1=self.n_CO_1, n_CO_2=self.n_CO_2, W_CO=self.CO.W,
                             Q1=self.Q1_(T, rho, Y_CH4, Y_O2), Q2=self.Q2_(T, rho, Y_CO, Y_O2, Y_CO2), degree=self.P)
        return rrCO

    def rrCO2_(self, T, rho, Y_CO, Y_O2, Y_CO2):
        """
        Compute the reaction rate of CO2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the reaction rate of CO2.
        """
        rrCO2 = do.Expression("W_CO2 * n_CO2_2 * Q2",
                              n_CO2_2=self.n_CO2_2, W_CO2=self.CO2.W,
                              Q2=self.Q2_(T, rho, Y_CO, Y_O2, Y_CO2), degree=self.P)
        return rrCO2

    def sO2dT_(self, T):
        """
        Compute the sensitivity of entropy of O2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of entropy of O2.
        """
        sO2dT = do.Expression("cp_O2 / T_",
                              cp_O2=self.janaf.janaf_cp_expr(self.O2, T), T_=T, degree=self.P)
        return sO2dT

    def sCOdT_(self, T):
        """
        Compute the sensitivity of entropy of CO with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of entropy of CO.
        """
        sCOdT = do.Expression("cp_CO / T_",
                              cp_CO=self.janaf.janaf_cp_expr(self.CO, T), T_=T, degree=self.P)
        return sCOdT

    def sCO2dT_(self, T):
        """
        Compute the sensitivity of entropy of CO2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of entropy of CO2.
        """
        sCO2dT = do.Expression("cp_CO2 / T_",
                               cp_CO2=self.janaf.janaf_cp_expr(self.CO2, T), T_=T, degree=self.P)
        return sCO2dT

    def Kf1dT_(self, T):
        """
        Compute the sensitivity of Kf1 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of Kf1.
        """
        Kf1dT = do.Expression("Kf1 * Ta1 / pow(T_, 2)",
                              Kf1=self.Kf1_(T), Ta1=self.Ta1, T_=T, degree=self.P)
        return Kf1dT

    def Kp2dT_(self, T):
        """
        Compute the sensitivity of Kp2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of Kp2.
        """
        Kp2dT = do.Expression("- Kp2 / (R * pow(T_, 2)) * (n_O2_2 * ha_O2 + n_CO_2 * ha_CO - n_CO2_2 * ha_CO2)",
                              Kp2=self.Kp2_(T), R=self.R, T_=T,
                              n_O2_2=self.n_O2_2, ha_O2=self.janaf.janaf_ha_expr(self.O2, T),
                              n_CO_2=self.n_CO_2, ha_CO=self.janaf.janaf_ha_expr(self.CO, T),
                              n_CO2_2=self.n_CO2_2, ha_CO2=self.janaf.janaf_ha_expr(self.CO2, T),
                              degree=self.P)
        return Kp2dT

    def Kc2dT_(self, T):
        """
        Compute the sensitivity of Kc2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of Kc2.
        """
        Kc2dT = do.Expression(
            "(Kp2dT - Kp2/T_ * (n_CO2_2 - n_O2_2 - n_CO_2)) * pow(pa / (R * T_), n_CO2_2-n_O2_2-n_CO_2)",
            Kp2dT=self.Kp2dT_(T), Kp2=self.Kp2_(T), pa=self.pa, R=self.R, T_=T,
            n_O2_2=self.n_O2_2, n_CO_2=self.n_CO_2, n_CO2_2=self.n_CO2_2,
            degree=self.P)
        return Kc2dT

    def Kf2dT_(self, T):
        """
        Compute the sensitivity of Kf2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of Kf2.
        """
        Kf2dT = do.Expression("Kf2 * Ta2 / pow(T_, 2)",
                              Kf2=self.Kf2_(T), Ta2=self.Ta2, T_=T, degree=self.P)
        return Kf2dT

    def Kr2dT_(self, T):
        """
        Compute the sensitivity of Kr2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.

        Returns
        -------
        Expression
            Expression for the sensitivity of Kr2.
        """
        Kr2dT = do.Expression("(Kf2dT * Kc2 - Kc2dT * Kf2) / pow(Kc2, 2)",
                              Kf2dT=self.Kf2dT_(T), Kc2dT=self.Kc2dT_(T), Kc2=self.Kc2_(T), Kf2=self.Kf2_(T),
                              degree=self.P)
        return Kr2dT

    def Q1dT_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the sensitivity of Q1 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q1.
        """
        Q1dT = do.Expression(
            "Kf1dT * pow(rho_, nu_CH4_1 + nu_O2_1) * pow(Y_CH4_ / W_CH4, nu_CH4_1) * pow(Y_O2_ / W_O2, nu_O2_1)",
            Kf1dT=self.Kf1dT_(T), rho_=rho, Y_CH4_=Y_CH4, Y_O2_=Y_O2,
            W_CH4=self.CH4.W, W_O2=self.O2.W,
            nu_CH4_1=self.nu_CH4_1, nu_O2_1=self.nu_O2_1, degree=self.P)
        return Q1dT

    def Q1drho_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the sensitivity of Q1 with respect to rho.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q1.
        """
        Q1drho = do.Expression("(nu_CH4_1 + nu_O2_1) / rho_ * Q1",
                               rho_=rho, nu_CH4_1=self.nu_CH4_1, nu_O2_1=self.nu_O2_1,
                               Q1=self.Q1_(T, rho, Y_CH4, Y_O2), degree=self.P)
        return Q1drho

    def Q1dYCH4_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the sensitivity of Q1 with respect to Y_CH4.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q1.
        """
        Y_CH4_thresh = do.Expression("Y_CH4_ < Y_CH4_lim_", Y_CH4_=Y_CH4, Y_CH4_lim_=self.YCH4_lim, degree=self.P)
        Y_CH4_limited = do.Expression("Y_CH4_ + Y_CH4_thresh_*(Y_CH4_lim_ - Y_CH4_)", Y_CH4_=Y_CH4,
                                      Y_CH4_thresh_=Y_CH4_thresh,
                                      Y_CH4_lim_=self.YCH4_lim, degree=self.P)
        Q1dYCH4 = do.Expression(
            "nu_CH4_1 / pow(W_CH4, nu_CH4_1) * Kf1 * pow(rho_, nu_CH4_1 + nu_O2_1) * pow(Y_CH4_, nu_CH4_1 - 1) * pow(Y_O2_ / W_O2, nu_O2_1)",
            rho_=rho, Y_CH4_=Y_CH4_limited, Y_O2_=Y_O2,
            Kf1=self.Kf1_(T), nu_CH4_1=self.nu_CH4_1, nu_O2_1=self.nu_O2_1,
            W_CH4=self.CH4.W, W_O2=self.O2.W, degree=self.P)
        return Q1dYCH4

    def Q1dYO2_(self, T, rho, Y_CH4, Y_O2):
        """
        Compute the sensitivity of Q1 with respect to Y_O2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CH4 : float
            Mass fraction of CH4.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q1.
        """
        Q1dYO2 = do.Expression("nu_O2_1 / Y_O2_ * Q1",
                               Y_O2_=Y_O2, nu_O2_1=self.nu_O2_1, Q1=self.Q1_(T, rho, Y_CH4, Y_O2), degree=self.P)
        return Q1dYO2

    def Q2dT_(self, T, rho, Y_CO, Y_O2, Y_CO2):
        """
        Compute the sensitivity of Q2 with respect to T.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q2.
        """
        Q2dT = do.Expression("Kf2dT * pow(rho_, nu_CO_2 + nu_O2_2) * pow(Y_CO_ / W_CO, nu_CO_2) * pow(Y_O2_ / W_O2, nu_O2_2) - \
                             Kr2dT * pow(rho_, nu_CO2_2) * pow(Y_CO2_ / W_CO2, nu_CO2_2)",
                             Kf2dT=self.Kf2dT_(T), Kr2dT=self.Kr2dT_(T),
                             rho_=rho, Y_CO_=Y_CO, Y_O2_=Y_O2, Y_CO2_=Y_CO2,
                             W_CO=self.CO.W, W_O2=self.O2.W, W_CO2=self.CO2.W,
                             nu_CO_2=self.nu_CO_2, nu_O2_2=self.nu_O2_2, nu_CO2_2=self.nu_CO2_2, degree=self.P)
        return Q2dT

    def Q2drho_(self, T, rho, Y_CO, Y_O2, Y_CO2):
        """
        Compute the sensitivity of Q2 with respect to rho.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q2.
        """
        Q2drho = do.Expression("(nu_CO_2 + nu_O2_2) / rho_ * Q2f - nu_CO2_2 / rho_ * Q2r",
                               Q2f=self.Q2f_(T, rho, Y_CO, Y_O2), Q2r=self.Q2r_(T, rho, Y_CO2), rho_=rho,
                               nu_CO_2=self.nu_CO_2, nu_O2_2=self.nu_O2_2, nu_CO2_2=self.nu_CO2_2, degree=self.P)
        return Q2drho

    def Q2dYCO_(self, T, rho, Y_O2):
        """
        Compute the sensitivity of Q2 with respect to Y_CO.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q2.
        """
        Q2dYCO = do.Expression(
            "nu_CO_2 / pow(W_CO, nu_CO_2) * Kf2 * pow(rho_, nu_CO_2 + nu_O2_2) * pow(Y_O2_ / W_O2, nu_O2_2)",
            Kf2=self.Kf2_(T), rho_=rho, Y_O2_=Y_O2, W_CO=self.CO.W, W_O2=self.O2.W,
            nu_CO_2=self.nu_CO_2, nu_O2_2=self.nu_O2_2, degree=self.P)
        return Q2dYCO

    def Q2dYO2_(self, T, rho, Y_CO, Y_O2):
        """
        Compute the sensitivity of Q2 with respect to Y_O2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q2.
        """
        Q2dYO2 = do.Expression("nu_O2_2 / Y_O2_ * Q2f",
                               Y_O2_=Y_O2, Q2f=self.Q2f_(T, rho, Y_CO, Y_O2), nu_O2_2=self.nu_O2_2, degree=self.P)
        return Q2dYO2

    def Q2dYCO2_(self, T, rho):
        """
        Compute the sensitivity of Q2 with respect to Y_CO2.

        Parameters
        ----------
        T : float
            Temperature.
        rho : float
            Density.
        Y_CO2 : float
            Mass fraction of CO2.

        Returns
        -------
        Expression
            Expression for the sensitivity of Q2.
        """
        Q2dYCO2 = do.Expression("- nu_CO2_2 / pow(W_CO2, nu_CO2_2) * Kr2 * pow(rho_, nu_CO2_2)",
                                rho_=rho, W_CO2=self.CO2.W, Kr2=self.Kr2_(T),
                                nu_CO2_2=self.nu_CO2_2, degree=self.P)
        return Q2dYCO2

    def Wmix(self, Y_CH4, Y_CO, Y_O2, Y_CO2, Y_H2O):
        """
        Compute the molecular weight of the mixture.

        Parameters
        ----------
        Y_CH4 : float
            Mass fraction of CH4.
        Y_CO : float
            Mass fraction of CO.
        Y_O2 : float
            Mass fraction of O2.
        Y_CO2 : float
            Mass fraction of CO2.
        Y_H2O : float
            Mass fraction of H2O.

        Returns
        -------
        float
            Molecular weight of the mixture.
        """
        Y_N2 = 1 - Y_CH4 - Y_CO - Y_O2 - Y_CO2 - Y_H2O
        return 1/(Y_CH4/self.CH4.W + Y_CO/self.CO.W + Y_O2/self.O2.W + Y_CO2/self.CO2.W + Y_H2O/self.H2O.W + Y_N2/self.N2.W)

    def add_source_to_weak_form(self, weakform, dQ_threshold=None):
        """
        Add a source term for the reaction and species to the weak form.

        Parameters
        ----------
        weakform : object
            The weak form to which the source term is added.
        dQ_threshold : float, optional
            Limits the reaction term to be applied only where the mean heat release rate is larger than this threshold.

        Returns
        -------
        Form
            Updated weak form with added source term.
        """
        if dQ_threshold:
            dQMean = weakform.dQMean
            reaction_active = do.project(do.Expression("dQ>dQmin", dQ=dQMean, dQmin=dQ_threshold, degree=weakform.order), 
                                         dQMean.function_space())
        else:
            reaction_active = do.Constant(1.0)

        if False:
            dQ1 = self.Q1dT * weakform.T + self.Q1drho * weakform.rho + self.Q1dYCH4 * weakform.YCH4 + self.Q1dYO2 * weakform.YO2
            dQ2 = self.Q2dT * weakform.T + self.Q2drho * weakform.rho + self.Q2dYCO * weakform.YCO + self.Q2dYCO2 * weakform.YCO2 + self.Q2dYO2 * weakform.YO2
            weakform.eq_energy_stiffness -= reaction_active * (- self.h0r_1 * dQ1 * weakform.v_ene - self.h0r_2 * dQ2 * weakform.v_ene) * weakform.dx
            weakform.eq_species_stiffness -= reaction_active * (- self.n_CH4_1 * self.CH4.W * dQ1 * weakform.v_YCH4) * weakform.dx
            weakform.eq_species_stiffness -= reaction_active * (- self.n_O2_1 * self.O2.W * dQ1 * weakform.v_YO2
                                                          - self.n_O2_2 * self.O2.W * dQ2 * weakform.v_YO2) * weakform.dx
            weakform.eq_species_stiffness -= reaction_active * (self.n_H2O_1 * self.H2O.W * dQ1 * weakform.v_YH2O) * weakform.dx
            weakform.eq_species_stiffness -= reaction_active * (self.n_CO_1 * self.CO.W * dQ1 * weakform.v_YCO
                                                          - self.n_CO_2 * self.CO.W * dQ2 * weakform.v_YCO ) * weakform.dx
            weakform.eq_species_stiffness -= reaction_active * (self.n_CO2_2 * self.CO2.W * dQ2 * weakform.v_YCO2) * weakform.dx

        else:
            dQ1R = self.Q1dT * weakform.TR + self.Q1drho * weakform.rhoR + self.Q1dYCH4 * weakform.YCH4R + self.Q1dYO2 * weakform.YO2R
            dQ2R = self.Q2dT * weakform.TR + self.Q2drho * weakform.rhoR + self.Q2dYCO * weakform.YCOR + self.Q2dYCO2 * weakform.YCO2R + self.Q2dYO2 * weakform.YO2R
            eq_stiffness = -reaction_active * (- self.h0r_1 * (dQ1R * weakform.v_eneR)
                                                         - self.h0r_2 * (dQ2R * weakform.v_eneR)) * weakform.dx
            eq_stiffness -= reaction_active * (- self.n_CH4_1 * self.CH4.W * (dQ1R * weakform.v_YCH4R)) * weakform.dx
            eq_stiffness -= reaction_active * (- self.n_O2_1 * self.O2.W * (dQ1R * weakform.v_YO2R)
                                                          - self.n_O2_2 * self.O2.W * (dQ2R * weakform.v_YO2R)) * weakform.dx
            eq_stiffness -= reaction_active * (self.n_H2O_1 * self.H2O.W * (dQ1R * weakform.v_YH2OR)) * weakform.dx
            eq_stiffness -= reaction_active * (self.n_CO_1 * self.CO.W * (dQ1R * weakform.v_YCOR)
                                                          - self.n_CO_2 * self.CO.W * (dQ2R * weakform.v_YCOR)) * weakform.dx
            eq_stiffness -= reaction_active * (self.n_CO2_2 * self.CO2.W * (dQ2R * weakform.v_YCO2R)) * weakform.dx

        return eq_stiffness
