import numpy as np
from dolfinx.fem import (
                        Expression,
                        )
class Janafopenfoam:

    def __init__(self, P):
        # order of basis functions, needed for expression degree
        self.P = P
        # universal gas constant
        self.R_univ = 8314.470066505450
        # standard temperature
        self.Tstd = 298.15
        # reference pressure
        self.pref = 1e5

    def janaf_cp(self, Specie, T):
        # Specieific isobaric heat-capacity cp [J/kmol K]
        # this is an implementation of the following OpenFOAM source-code:
        # https://github.com/OpenFOAM/OpenFOAM-2.3.x/blob/master/src/thermophysicalModels/Specieie/thermo/janaf/janafThermoI.H
        if (T < Specie.Tlow) or (T > Specie.Thigh):
            raise Exception("Temperature not in temperature range for JANAF polynomial.")

        if T < Specie.Tcommon:
            cp = self.R_univ * ((((Specie.lowCpCoeffs[4]*T + Specie.lowCpCoeffs[3])*T
                                  + Specie.lowCpCoeffs[2])*T + Specie.lowCpCoeffs[1])*T + Specie.lowCpCoeffs[0])
        else:
            cp = self.R_univ * ((((Specie.highCpCoeffs[4]*T + Specie.highCpCoeffs[3])*T
                                  + Specie.highCpCoeffs[2])*T + Specie.highCpCoeffs[1])*T + Specie.highCpCoeffs[0])
        return cp

    def janaf_cp_expr(self, Specie, T):
        # Specieific isobaric heat-capacity cp [J/kmol K]
        # computed by expression to be used on dolfin fields
        # this is an implementation of the following OpenFOAM source-code:
        # https://github.com/OpenFOAM/OpenFOAM-2.3.x/blob/master/src/thermophysicalModels/Specieie/thermo/janaf/janafThermoI.H

        # check where the high temperature coefficients are needed
        # returns 1 where T > Tcom, 0 elsewhere
        highCoeff = Expression("T_ > Tcom", T_=T, Tcom=Specie.Tcommon, degree=self.P)

        # compute cp like this: use the low coefficients by default, ad the difference between
        # high and low coefficients, where hC is equal to 1
        cp = do.Expression("R * ((((lc4*T_ + lc3)*T_ + lc2)*T_ + lc1)*T_ + lc0) + highCoeff_ * R * (((((hc4-lc4)*T_ + (hc3-lc3))*T_ + (hc2-lc2))*T_ + (hc1-lc1))*T_ + (hc0-lc0))",
                           R=self.R_univ, T_=T, highCoeff_=highCoeff,
                           lc0=Specie.lowCpCoeffs[0], lc1=Specie.lowCpCoeffs[1], lc2=Specie.lowCpCoeffs[2],
                           lc3=Specie.lowCpCoeffs[3], lc4=Specie.lowCpCoeffs[4],
                           hc0=Specie.highCpCoeffs[0], hc1=Specie.highCpCoeffs[1], hc2=Specie.highCpCoeffs[2],
                           hc3=Specie.highCpCoeffs[3], hc4=Specie.highCpCoeffs[4],
                           degree=self.P)
        return cp

    def janaf_hc(self, Specie):
        # chemical enthalpy hc, also called enthalpy of formation [J/kmol]
        hc = self.R_univ * (((((Specie.lowCpCoeffs[4]/5*self.Tstd + Specie.lowCpCoeffs[3]/4)*self.Tstd
                             + Specie.lowCpCoeffs[2]/3)*self.Tstd + Specie.lowCpCoeffs[1]/2)*self.Tstd
                             + Specie.lowCpCoeffs[0])*self.Tstd + Specie.lowCpCoeffs[5])
        return hc

    def janaf_ha(self, Specie, T):
        # absolute enthalpy ha [J/kmol]
        if (T < Specie.Tlow) or (T > Specie.Thigh):
            raise Exception("Temperature not in temperature range for JANAF polynomial.")

        if T < Specie.Tcommon:
            ha = self.R_univ * (((((Specie.lowCpCoeffs[4]/5*T + Specie.lowCpCoeffs[3]/4)*T + Specie.lowCpCoeffs[2]/3)*T
                                  + Specie.lowCpCoeffs[1]/2)*T + Specie.lowCpCoeffs[0])*T + Specie.lowCpCoeffs[5])
        else:
            ha = self.R_univ * (((((Specie.highCpCoeffs[4]/5*T + Specie.highCpCoeffs[3]/4)*T + Specie.highCpCoeffs[2]/3)*T
                                  + Specie.highCpCoeffs[1]/2)*T + Specie.highCpCoeffs[0])*T + Specie.highCpCoeffs[5])
        return ha

    def janaf_ha_expr(self, Specie, T):
        # absolute enthalpy ha [J/kmol]
        # computed by expression to be used on dolfin fields
        # check where the high temperature coefficients are needed
        # returns 1 where T > Tcom, 0 elsewhere
        highCoeff = Expression("T_ > Tcom", T_=T, Tcom=Specie.Tcommon, degree=self.P)

        # compute ha like this: use the low coefficients by default, ad the difference between
        # high and low coefficients, where hC is equal to 1
        ha = Expression("R * (((((lc4/5*T_ + lc3/4)*T_ + lc2/3)*T_ + lc1/2)*T_ + lc0)*T_ + lc5) + highCoeff_ * R * ((((((hc4-lc4)/5*T_ + (hc3-lc3)/4)*T_ + (hc2-lc2)/3)*T_ + (hc1-lc1)/2)*T_ + (hc0-lc0))*T_ + (hc5-lc5))",
                           R=self.R_univ, T_=T, highCoeff_=highCoeff,
                           lc0=Specie.lowCpCoeffs[0], lc1=Specie.lowCpCoeffs[1], lc2=Specie.lowCpCoeffs[2], lc3=Specie.lowCpCoeffs[3],
                           lc4=Specie.lowCpCoeffs[4], lc5=Specie.lowCpCoeffs[5],
                           hc0=Specie.highCpCoeffs[0], hc1=Specie.highCpCoeffs[1], hc2=Specie.highCpCoeffs[2], hc3=Specie.highCpCoeffs[3],
                           hc4=Specie.highCpCoeffs[4], hc5=Specie.highCpCoeffs[5],
                           degree=self.P)

        return ha

    def janaf_hs(self, Specie, T):
        # sensible enthalpy hs [J/kmol]
        if (T < Specie.Tlow) or (T > Specie.Thigh):
            raise Exception("Temperature not in temperature range for JANAF polynomial.")
        hs = self.janaf_ha(Specie, T) - self.janaf_hc(Specie)
        return hs

    def janaf_hs_expr(self, Specie, T):
        # sensible enthalpy hs [J/kmol]
        # computed by expression to be used on dolfin fields
        hs = Expression("ha - hc", ha=self.janaf_ha_expr(Specie, T), hc=self.janaf_hc(Specie), degree=self.P)
        return hs

    def janaf_Hc(self, Specie):
        # chemical enthalpy Hc [J/kg]
        Hc = self.janaf_hc(Specie) / Specie.W
        return Hc

    def janaf_Ha(self, Specie, T):
        # absolute enthalpy Ha [J/kg]
        Ha = self.janaf_ha(Specie, T) / Specie.W
        return Ha

    def janaf_Ha_expr(self, Specie, T):
        # absolute enthalpy Ha [J/kg]
        # computed by expression to be used on dolfin fields
        Ha = Expression("ha / W", ha=self.janaf_ha_expr(Specie, T), W=Specie.W, degree=self.P)
        return Ha

    def janaf_Hs(self, Specie, T):
        # sensible enthalpy Hs [J/kg]
        Hs = self.janaf_hs(Specie, T) / Specie.W
        return Hs

    def janaf_Hs_expr(self, Specie, T):
        # sensible enthalpy Hs [J/kg]
        # computed by expression to be used on dolfin fields
        Hs = Expression("hs / W", hs=self.janaf_hs_expr(Specie, T), W=Specie.W, degree=self.P)
        return Hs

    def janaf_s0(self, Specie):
        # standard entropy s0 [J/kmol K]
        # this is an implementation of the following OpenFOAM source-code:
        # https://github.com/OpenFOAM/OpenFOAM-2.3.x/blob/master/src/thermophysicalModels/specie/thermo/janaf/janafThermoI.H
        s0 = self.R_univ * ((((Specie.lowCpCoeffs[4]/4*self.Tstd + Specie.lowCpCoeffs[3]/3)*self.Tstd
                              + Specie.lowCpCoeffs[2]/2)*self.Tstd + Specie.lowCpCoeffs[1])*self.Tstd
                            + Specie.lowCpCoeffs[0]*np.log(self.Tstd) + Specie.lowCpCoeffs[6])
        return s0

    def janaf_s(self, Specie, T):
        # entropy s [J/kmol K]
        if (T < Specie.Tlow) or (T > Specie.Thigh):
            raise Exception("Temperature not in temperature range for JANAF polynomial.")

        if T < Specie.Tcommon:
            s = self.R_univ * ((((Specie.lowCpCoeffs[4]/4*T + Specie.lowCpCoeffs[3]/3)*T
                                 + Specie.lowCpCoeffs[2]/2)*T + Specie.lowCpCoeffs[1])*T
                               + Specie.lowCpCoeffs[0]*np.log(T) + Specie.lowCpCoeffs[6])
        else:
            s = self.R_univ * ((((Specie.highCpCoeffs[4]/4*T + Specie.highCpCoeffs[3]/3)*T
                                 + Specie.highCpCoeffs[2]/2)*T + Specie.highCpCoeffs[1])*T
                               + Specie.highCpCoeffs[0]*np.log(T) + Specie.highCpCoeffs[6])
        return s

    def janaf_s_expr(self, Specie, T):
        # entropy s [J/kmol K]
        # computed by expression to be used on dolfin fields
        # check where the high temperature coefficients are needed
        # returns 1 where T > Tcom, 0 elsewhere
        highCoeff = Expression("T_ > Tcom", T_=T, Tcom=Specie.Tcommon, degree=self.P)

        # compute ha like this: use the low coefficients by default, ad the difference between
        # high and low coefficients, where hC is equal to 1
        # use the change of base for logarithms to compute the natural logarithm with log10 because I couldn't find out how ln works with Expressions ...
        s = Expression("R * ((((lc4/4*T_ + lc3/3)*T_ + lc2/2)*T_ + lc1)*T_ + lc0*log10(T_)/log10(e) + lc6) + highCoeff_ * R * (((((hc4-lc4)/4*T_ + (hc3-lc3)/3)*T_ + (hc2-lc2)/2)*T_ + (hc1-lc1))*T_ + (hc0-lc0)*log10(T_)/log10(e) + (hc6-lc6))",
                          R=self.R_univ, T_=T, highCoeff_=highCoeff, e=np.e,
                          lc0=Specie.lowCpCoeffs[0], lc1=Specie.lowCpCoeffs[1], lc2=Specie.lowCpCoeffs[2], lc3=Specie.lowCpCoeffs[3],
                          lc4=Specie.lowCpCoeffs[4], lc6=Specie.lowCpCoeffs[6],
                          hc0=Specie.highCpCoeffs[0], hc1=Specie.highCpCoeffs[1], hc2=Specie.highCpCoeffs[2], hc3=Specie.highCpCoeffs[3],
                          hc4=Specie.highCpCoeffs[4], hc6=Specie.highCpCoeffs[6],
                          degree=self.P)
        return s

    def janaf_S0(self, Specie):
        # standard Entropy S0 [J/kg K]
        S0 = self.janaf_s0(Specie) / Specie.W
        return S0

    def janaf_S(self, Specie, T):
        # Entropy S [J/kg K]
        S = self.janaf_s(Specie, T) / Specie.W
        return S

    def janaf_S_expr(self, Specie, T):
        # standard Entropy S [J/kg K]
        # computed by expression to be used on dolfin fields
        S = Expression("s / W", s=self.janaf_s_expr(Specie, T), W=Specie.W, degree=self.P)
        return S

    def janaf_S_in_mix_expr(self, Specie, T, p, Y_spec, Wmix):
        # Entropy S of one species in a mixture [J/kg K]
        # computed by expression to be used on dolfin fields

        # account for possible log(0) if Y_spec == 0 by thresholding
        thresh = 1e-6
        Y_spec_thresh = Expression("Y_spec < thresh", Y_spec=Y_spec, thresh=thresh, degree=self.P)
        Y_spec_limited = Expression("Y_spec + Y_spec_thresh*(thresh - Y_spec)", Y_spec=Y_spec, Y_spec_thresh=Y_spec_thresh, thresh=thresh, degree=self.P)

        S_in_mix = Expression("S - R * log10(Y_spec*Wmix/W*p/pref)/log10(e)",
                                 S=self.janaf_S_expr(Specie, T), R=Specie.R, Y_spec=Y_spec_limited,
                                 Wmix=Wmix, p=p, pref=self.pref, W=Specie.W, e=np.e, degree=self.P)
        return S_in_mix

    def janaf_g(self, Specie, T):
        # gibbs free energy g [J/kmol]
        g = self.janaf_ha(Specie, T) - T * self.janaf_s(Specie, T)
        return g

    def janaf_gs(self, Specie, T):
        # sensible gibbs free energy g [J/kmol]
        gs = self.janaf_hs(Specie, T) - T * self.janaf_s(Specie, T)
        return gs

    def janaf_g_expr(self, Specie, T):
        # gibbs free energy g [J/kmol]
        # computed by expression to be used on dolfin fields
        g = Expression("ha - T_ * s", ha=self.janaf_ha_expr(Specie, T), T_=T, s=self.janaf_s_expr(Specie, T), degree=self.P)
        return g

    def janaf_gs_expr(self, Specie, T):
        # sensible gibbs free energy g [J/kmol]
        # computed by expression to be used on dolfin fields
        gs = Expression("hs - T_ * s", hs=self.janaf_hs_expr(Specie, T), T_=T, s=self.janaf_s_expr(Specie, T), degree=self.P)
        return gs

    def janaf_G(self, Specie, T):
        # gibbs free energy G [J/kg]
        G = self.janaf_g(Specie, T) / Specie.W
        return G

    def janaf_Gs(self, Specie, T):
        # sensible gibbs free energy G [J/kg]
        Gs = self.janaf_gs(Specie, T) / Specie.W
        return Gs

    def janaf_G_expr(self, Specie, T):
        # gibbs free energy G [J/kg]
        # computed by expression to be used on dolfin fields
        G = Expression("g / W", gs=self.janaf_g_expr(Specie, T), W=Specie.W, degree=self.P)
        return G

    def janaf_Gs_expr(self, Specie, T):
        # sensible gibbs free energy G [J/kg]
        # computed by expression to be used on dolfin fields
        Gs = Expression("gs / W", gs=self.janaf_gs_expr(Specie, T), W=Specie.W, degree=self.P)
        return Gs
