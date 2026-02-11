#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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



from Tabs.JANAFTabulate1971 import getEnthalpyPlusFormationJANAF1971,getEntropyJANAF1971
from dolfinx.fem import (
                        Expression,
                        Function,
                        )
from ufl import (
                dx,
)
import numpy as np
import bisect as bs
class TwoStepReaction():
    """
    Represents a two-step methane-air combustion reaction mechanism.

    This class models two-step chemical reactions, adapted for BFER or 2S_CH4_CM2 models, and provides methods for reading reaction data, computing mean fields, and evaluating reaction rates.
    See CERFACS for model reference https://www.cerfacs.fr/cantera/mechanisms/meth.php
    
    **Initialize the TwoStepReaction object**

    Parameters
    ----------
    ModelName : str
        Name of the reaction model to use.

    Attributes
    ----------
    mixtureDirectory : str
        Directory name for mixture data.
    speciesDirectory : str
        Directory name for species data.
    p0 : float
        Reference pressure.
    R : float
        Universal gas constant.
    n_CH4_1 : float or None
        Stoichiometric prefactor for CH4 in reaction 1.
    n_O2_1 : float or None
        Stoichiometric prefactor for O2 in reaction 1.
    n_CO_1 : float or None
        Stoichiometric prefactor for CO in reaction 1.
    n_CO_2 : float or None
        Stoichiometric prefactor for CO in reaction 2.
    n_O2_2 : float or None
        Stoichiometric prefactor for O2 in reaction 2.
    n_CO2_2 : float or None
        Stoichiometric prefactor for CO2 in reaction 2.
    nu_CH4_1 : float or None
        Kinematic exponent for CH4 in reaction 1.
    nu_O2_1 : float or None
        Kinematic exponent for O2 in reaction 1.
    nu_CO_2 : float or None
        Kinematic exponent for CO in reaction 2.
    nu_O2_2 : float or None
        Kinematic exponent for O2 in reaction 2.
    nu_CO2_2 : float or None
        Kinematic exponent for CO2 in reaction 2.
    A1 : float or None
        Pre-exponential factor for reaction 1.
    Ta1 : float or None
        Activation temperature for reaction 1.
    beta1 : float or None
        Temperature exponent for reaction 1.
    h1 : float or None
        Heat of reaction for reaction 1.
    A2 : float or None
        Pre-exponential factor for reaction 2.
    Ta2 : float or None
        Activation temperature for reaction 2.
    beta2 : float or None
        Temperature exponent for reaction 2.
    h2 : float or None
        Heat of reaction for reaction 2.
    WO2 : float or None
        Molar mass of O2.
    WH2O : float or None
        Molar mass of H2O.
    WCH4 : float or None
        Molar mass of CH4.
    WCO : float or None
        Molar mass of CO.
    WCO2 : float or None
        Molar mass of CO2.
    WN2 : float or None
        Molar mass of N2.
    Q1 : object or None
        Reaction progress rate for reaction 1.
    Q2f : object or None
        Forward reaction progress rate for reaction 2.
    Q2r : object or None
        Reverse reaction progress rate for reaction 2.
    Q2 : object or None
        Total reaction progress rate for reaction 2.
    K2f : object or None
        Forward reaction rate constant for reaction 2.
    K2r : object or None
        Reverse reaction rate constant for reaction 2.
    lnexpInEqui_h : object or None
        Linearized enthalpy term for equilibrium constant in reaction 2.
    expInEqui : object or None
        Exponential term for equilibrium constant in reaction 2.
    Equi : object or None
        Equilibrium constant for reaction 2.
    i_rho : int or None
        Index for density in solution variables.
    i_CH4 : int or None
        Index for CH4 in solution variables.
    i_O2 : int or None
        Index for O2 in solution variables.
    i_CO2 : int or None
        Index for CO2 in solution variables.
    i_CO : int or None
        Index for CO in solution variables.
    T : object or None
        Temperature field.
    Y_CH4_limited : object or None
        Limited CH4 mass fraction.
    epsilon : float
        Correction for base flow CH4.
    reactionName : str
        Name of the reaction.
    """

    def __init__(self,
    ModelName,
    ):
        """
        Initializes the TwoStepReaction instance.

        Parameters
        ----------
        ModelName : str
            Name of the reaction model to use.
        """

        self.mixtureDirectory = 'Mixture' #to be put in param
        self.speciesDirectory = 'Species' #to be put in param



        self.p0=101300
        self.R=8.31446261815324

        # stochiometric prefactors
        self.n_CH4_1 = None
        self.n_O2_1 = None
        self.n_CO_1 = None
        #self.n_H2O_1 = 2

        self.n_CO_2 = None
        self.n_O2_2 = None
        self.n_CO2_2 = None


        # kinematic exponents:
        self.nu_CH4_1 = None
        self.nu_O2_1 = None

        self.nu_CO_2 = None
        self.nu_O2_2 = None
        self.nu_CO2_2 = None

        # reaction constants
        self.A1 = None
        self.Ta1 = None
        self.beta1 = None

        self.h1 = None


        self.A2 = None
        self.Ta2 = None
        self.beta2 = None

        self.h2 = None

        # molar mass
        self.WO2 = None
        self.WH2O = None
        self.WCH4 = None
        self.WCO = None
        self.WCO2 = None
        self.N2 = None


        # reaction progress rates
        self.Q1 = None #total reaction 1

        self.Q2f = None #forwared reaction in reaction 2
        self.Q2r = None #inverse reaction in reaction 2
        self.Q2 = None #total reaction 2

        # reaction rate (only as function of temperature)
        self.K2f=None
        self.K2r=None

        # variables to be linearized step-by-step in equilibrium contant in reaction 2
        self.lnexpInEqui_h=None
        self.expInEqui=None
        self.Equi=None


        self.i_rho=None
        self.i_CH4=None
        self.i_O2=None
        self.i_CO2=None
        self.i_CO=None
        self.T=None

#       #PEA fitting functions (homogeneously equal to 1 evaluated for lean premixed case.)
        #if required, formula and parameters given at the end of script)
#       self.phi=None
#       self.f1=None
#       self.f2=None


        self.Y_CH4_limited=None

        self.epsilon=0.0002 #correction for base flow CH4


        self.reactionName='TwoStep'
        print('Initializing reaction '+self.reactionName+' '+ModelName)
        self.read_reaction_dict(ModelName)


    def read_reaction_dict(self,
    ModelName,
    ):
        """
        Reads reaction parameters and species data from files and updates class attributes.

        Parameters
        ----------
        ModelName : str
            Name of the reaction model to use.

        Notes
        -----
        This method sets stoichiometric prefactors, kinematic exponents, reaction constants, and molar masses based on the provided model name.
        """
        fileMixture = open(
        self.mixtureDirectory,
        'r',
        )
        mixtureDictDict = eval(fileMixture.read())
        fileMixture.close()
        MD = mixtureDictDict[ModelName]
        #print(MD)
        # stochiometric prefactors
        self.n_CH4_1 = MD['reac_st_CH4_1']
        self.n_O2_1 = MD['reac_st_O2_1']
        self.n_CO_1 = MD['reac_st_CO_1']

        self.n_CO_2 = MD['reac_st_CO_2']
        self.n_O2_2 = MD['reac_st_O2_2']
        self.n_CO2_2 = MD['reac_st_CO2_2']


        # kinematic exponents:
        self.nu_CH4_1 = MD['reac_nu_CH4_1']
        self.nu_O2_1 = MD['reac_nu_O2_1']

        self.nu_CO_2 = MD['reac_nu_CO_2']
        self.nu_O2_2 = MD['reac_nu_O2_2']
        self.nu_CO2_2 = MD['reac_nu_CO2_2']

        # reaction constants
        self.A1 = MD['reac_preexp_1']
        self.Ta1 = MD['reac_act_tem_1']
        self.beta1 = MD['reac_exp_tem_1']

        self.h1 = MD['reac_h0_1']

        self.A2 = MD['reac_preexp_2']
        self.Ta2 = MD['reac_act_tem_2']
        self.beta2 = MD['reac_exp_tem_2']

        self.h2 = MD['reac_h0_2']

        fileSpecies = open(
        self.speciesDirectory,
        'r',
        )
        speciesDictDict = eval(fileSpecies.read())
        fileSpecies.close()
        self.WCH4 = speciesDictDict['CH4']['mol_weight']
        self.WO2 = speciesDictDict['O2']['mol_weight']
        self.WCO2 = speciesDictDict['CO2']['mol_weight']
        self.WCO = speciesDictDict['CO']['mol_weight']
        self.WH2O = speciesDictDict['H2O']['mol_weight']
        self.WN2 = speciesDictDict['N2']['mol_weight']

    def compute_mean_field(self,
    MF,
    ele,
    ):
        """
        Computes the mean field reaction rates and equilibrium constants.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Returns
        -------
        Q1 : object
            Mean field for reaction 1.
        Q2f : object
            Forward mean field for reaction 2.
        Q2r : object
            Reverse mean field for reaction 2.

        Notes
        -----
        This method also computes equilibrium constants and related fields for reaction 2.
        """
        self.Y_CH4_limited = Expression(
        "Y_CH4_ + epsilon_",
        Y_CH4_=MF['CH4'],
        epsilon_=self.epsilon,
        degree=2,
        )

        ele.interpolate(Expression(
        "A*exp(-Ta/T_) * pow(T_, beta)* pow(rho_, nu_CH4 + nu_O2) * pow(Y_CH4_ / W_CH4, nu_CH4) * pow(Y_O2_ / W_O2, nu_O2)",
        A=self.A1,
        T_=MF['T'],
        Ta=self.Ta1,
        beta=self.beta1,
        rho_=MF['rho'],
        nu_CH4=self.nu_CH4_1,
        nu_O2=self.nu_O2_1,
        Y_CH4_=self.Y_CH4_limited,
        W_CH4=self.WCH4,
        Y_O2_=MF['O2'],
        W_O2=self.WO2,
        degree=2,
        ))

        self.Q1 = ele
        self.K2f=project(
        Expression(
        "A*exp(-Ta/T_) * pow(T_, beta)",
        A=self.A2,
        T_=MF['T'],
        Ta=self.Ta2,
        beta=self.beta2,
        degree=2,
        ),
        ele,
        )
        self.Q2f=project(
        Expression(
        "K2f_*pow(rho_, nu_CO + nu_O2) * pow(Y_CO_ / W_CO, nu_CO) * pow(Y_O2_ / W_O2, nu_O2)",
        K2f_=self.K2f,
        rho_=MF['rho'],
        nu_CO=self.nu_CO_2,
        nu_O2=self.nu_O2_2,
        Y_CO_=MF['CO'],
        W_CO=self.WCO,
        Y_O2_=MF['O2'],
        W_O2=self.WO2,
        degree=2,
        ),
        ele,
        )

        #choose the base flow of temperature imported from AVBP. This can reproduce more precisely the base flow reaction rates
        self.T=MF['T']
        Tfield=MF['T'].vector[:]


        lnexpInEqui_h = ((-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'CO2',
        Tfield,
        )))/self.R/Tfield)*self.n_CO2_2\
                -(-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'CO',
        Tfield,
        )))/self.R/Tfield)*self.n_CO_2\
                -(-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'O2',
        Tfield,
        )))/self.R/Tfield)*self.n_O2_2)
        lnexpInEqui = (np.array(list(getEntropyJANAF1971(
        'CO2',
        Tfield,
        )))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'CO2',
        Tfield,
        )))/self.R/Tfield)*self.n_CO2_2\
                -(np.array(list(getEntropyJANAF1971(
        'CO',
        Tfield,
        )))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'CO',
        Tfield,
        )))/self.R/Tfield)*self.n_CO_2\
                -(np.array(list(getEntropyJANAF1971(
        'O2',
        Tfield,
        )))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971(
        'O2',
        Tfield,
        )))/self.R/Tfield)*self.n_O2_2

        nu_j=self.nu_CO2_2-self.nu_CO_2-self.nu_O2_2
        self.lnexpInEqui_h=Function(ele)
        self.lnexpInEqui_h.vector[:]=lnexpInEqui_h
        self.expInEqui=Function(ele)
        self.expInEqui.vector[:]=np.exp(lnexpInEqui)
        self.Equi=Function(ele)
        self.Equi.vector[:]=(self.p0/self.R/Tfield)**nu_j*np.exp(lnexpInEqui)

        self.K2r=project(
        Expression(
        "k2f_/Equi_",
        k2f_=self.K2f,
        Equi_=self.Equi,
        degree=2,
        ),
        ele,
        )
        self.Q2r=project(
        Expression(
        "K2r_*pow(rho_, nu_CO2) * pow(Y_CO2_ / W_CO2, nu_CO2) ",
        K2r_=self.K2r,
        rho_=MF['rho'],
        nu_CO2=self.nu_CO2_2,
        Y_CO2_=MF['CO2'],
        W_CO2=self.WCO2,
        degree=2,
        ),
        ele,
        )



        return self.Q1, self.Q2f, self.Q2r


    def add_reaction(self,
    MF,
    testf,
    trialf,
    solutionList,
    ele,
    ):
        """
        Adds the reaction terms to the weak form for both reaction steps.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        testf : list
            List of test functions.
        trialf : dict
            Dictionary of trial functions.
        solutionList : list of str
            List of solution variable names.
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Returns
        -------
        reaction_term : ufl.Form
            Weak form expression representing the reaction contributions.

        Notes
        -----
        The method computes the indices for relevant variables and constructs the weak form using the reaction rates and heat release for both steps.
        """
        self.i_rho=solutionList.index('rho')
        self.i_CH4=solutionList.index('CH4')
        self.i_O2=solutionList.index('O2')
        self.i_CO2=solutionList.index('CO2')
        self.i_CO=solutionList.index('CO')

        #dQ1
        dQ1=self.d_q1_(
        MF,
        trialf,
        )
        form = -testf[self.i_rho]*dQ1*self.Q1*self.h1*dx\
            -testf[self.i_CH4]*dQ1*self.Q1*self.n_CH4_1*self.WCH4*dx\
            -testf[self.i_O2]*dQ1*self.Q1*self.n_O2_1*self.WO2*dx\
            +testf[self.i_CO]*dQ1*self.Q1*self.n_CO_1*self.WCO*dx


#       #dQ2f
        dQ2f=self.d_q2f_(
        MF,
        trialf,
        )
        dQ2r=self.d_q2r_(
        MF,
        trialf,
        )
        dQ2=dQ2f*self.Q2f-dQ2r*self.Q2r

        form += \
            -testf[self.i_O2]*dQ2*self.n_O2_2*self.WO2*dx\
            -testf[self.i_CO]*dQ2*self.n_CO_2*self.WCO*dx\
            +testf[self.i_CO2]*dQ2*self.n_CO2_2*self.WCO2*dx\
            -testf[self.i_rho]*dQ2*self.h2*dx
        return form

    def d_q1_(self,
    MF,
    trialf,
    ):
        """
        Computes the fluctuation of the reaction rate for the first reaction step.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        trialf : dict
            Dictionary of trial functions.

        Returns
        -------
        dQ1 : float
            Fluctuation of the reaction rate for reaction 1.
        """
        dT = -trialf['rho']/MF['rho']*MF['T']
        return ((self.nu_O2_1+self.nu_CH4_1)*trialf['rho']/MF['rho']\
                +self.beta1*dT/MF['T']\
                +self.Ta1*dT/MF['T']/MF['T']\
                +self.nu_O2_1*trialf['O2']/MF['O2']\
                +self.nu_CH4_1*trialf['CH4']/(MF['CH4']+self.epsilon))#((MF['O2']-0.0447)/1.78e-1*4.45e-2))#/MF['CH4'])#self.Y_CH4_limited)


    def d_q2f_(self,
    MF,
    trialf,
    ):
        """
        Computes the fluctuation of the forward reaction rate for the second reaction step.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        trialf : dict
            Dictionary of trial functions.

        Returns
        -------
        dQ2f : float
            Fluctuation of the forward reaction rate for reaction 2.
        """
        dT = -trialf['rho']/MF['rho']*MF['T']
        return ((self.nu_O2_2+self.nu_CO_2)*trialf['rho']/MF['rho']\
                +self.beta2*dT/MF['T']\
                +self.Ta2*dT/MF['T']/MF['T']\
                +self.nu_CO_2*trialf['CO']/(MF['CO'])\
                +self.nu_O2_2*trialf['O2']/MF['O2'])

    def d_q2r_(self,
    MF,
    trialf,
    ):
        """
        Computes the fluctuation of the reverse reaction rate for the second reaction step.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        trialf : dict
            Dictionary of trial functions.

        Returns
        -------
        dQ2r : float
            Fluctuation of the reverse reaction rate for reaction 2.
        """
        dK2rdT=self.d_k2rd_t_(self.T)
        dT = -trialf['rho']/MF['rho']*MF['T']
        return (\
                +self.nu_CO2_2*trialf['CO2']/(MF['CO2'])\
                +self.nu_CO2_2*trialf['rho']/MF['rho']\
                +dT*dK2rdT\
                )
#       #if we neglect dK2rdT
#       return (\
#               +self.nu_CO2_2*trialf['CO2']/(MF['CO2'])\
#               +self.nu_CO2_2*trialf['rho']/MF['rho']\
#               +dT*dK2rdT\
#               )

    def d_k2fd_t_(self,
    T,
    ):
        """
        Computes the temperature derivative of the forward reaction rate constant for reaction 2.

        Parameters
        ----------
        T : object
            Temperature field.

        Returns
        -------
        dK2fdT : object
            Expression for the temperature derivative of the forward rate constant.
        """
        return Expression(
        "(Ta2 / pow(T_, 2)+beta2/T_)",
        Ta2=self.Ta2,
        T_=T,
        beta2=self.beta2,
        degree=2,
        )

    def d_k2rd_t_(self,
    T,
    ):
        """
        Computes the temperature derivative of the reverse reaction rate constant for reaction 2.

        Parameters
        ----------
        T : object
            Temperature field.

        Returns
        -------
        dK2rdT : object
            Expression for the temperature derivative of the reverse rate constant.
        """

        return Expression(
        "Kf2dT_ -KcdT_/Kc",
        kf2d_t_=self.d_k2fd_t_(T),
        Kc=self.Equi,
        K2f_=self.K2f,
        K2r_=self.K2r,
        KcdT_=self.d_equid_t_(self.T),
        degree=2,
        )
#       #if we neglect KcdT
#       return Expression("Kf2dT_ ",Kf2dT_=self.dK2fdT_(T), Kc=self.Equi, K2f_=self.K2f, K2r_=self.K2r,\
#                   KcdT_=self.dEquidT_(self.T),  degree=2)

    def d_equid_t_(self,
    T,
    ):
        """
        Computes the temperature derivative of the equilibrium constant for reaction 2.

        Parameters
        ----------
        T : object
            Temperature field.

        Returns
        -------
        dEquidT : object
            Expression for the temperature derivative of the equilibrium constant.
        """
        return Expression(
        "(Kp2dT - Kp2/T_ * (n_CO2_2 - n_O2_2 - n_CO_2)) * pow(pa / (R * T_), n_CO2_2-n_O2_2-n_CO_2)",
        Kp2dT=self.dexp_in_equid_t_(T),
        Kp2=self.expInEqui,
        pa=self.p0,
        R=self.R,
        T_=T,
        n_O2_2=self.n_O2_2,
        n_CO_2=self.n_CO_2,
        n_CO2_2=self.n_CO2_2,
        degree=2,
        )
#       #if we neglect Kp2dT
#       return Expression(
#       "(- Kp2/T_ * (n_CO2_2 - n_O2_2 - n_CO_2)) * pow(pa / (R * T_), n_CO2_2-n_O2_2-n_CO_2)",
#        Kp2=self.expInEqui, pa=self.p0, R=self.R, T_=T,
#       n_O2_2=self.n_O2_2, n_CO_2=self.n_CO_2, n_CO2_2=self.n_CO2_2,
#       degree=2)

    def dexp_in_equid_t_(self,
    T,
    ):
        """
        Computes the temperature derivative of the exponential part of the equilibrium constant for reaction 2.

        Parameters
        ----------
        T : object
            Temperature field.

        Returns
        -------
        dexpInEquidT : object
            Expression for the temperature derivative of the exponential part of the equilibrium constant.
        """
        return Expression(
        "- Kp2 /  T_ * H_",
        Kp2=self.expInEqui,
        R=self.R,
        T_=T,
        H_=self.lnexpInEqui_h,
        degree=2,
        )


    def post_heat_relaese(self,
    MF,
    prho,
    pCH4,
    pO2,
    pCO,
    pCO2,
    ele,
    ):
        """
        Evaluates the fluctuation of heat release in post-processing.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        prho : float
            Perturbation in density.
        pCH4 : float
            Perturbation in CH4 mass fraction.
        pO2 : float
            Perturbation in O2 mass fraction.
        pCO : float
            Perturbation in CO mass fraction.
        pCO2 : float
            Perturbation in CO2 mass fraction.
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Returns
        -------
        heat_release_field : object
            Projected fluctuation heat release field.
        """
        #evalute fluctuation HeatRelease in post-processing
        #dQ1
        dQ1=self.postd_q1_(
        MF,
        prho,
        pO2,
        pCH4,
        ele,
        )
        form = dQ1*self.Q1*self.h1

        #dQ2
        dQ2f=self.postd_q2f_(
        MF,
        prho,
        pCO,
        pO2,
        )
        dQ2r=self.postd_q2r_(
        MF,
        prho,
        pCO2,
        )
        form += (dQ2f*self.Q2f-dQ2r*self.Q2r)*self.h2

        return project(
        form,
        ele,
        )

    def postd_q1_(self,
    MF,
    prho,
    pO2,
    pCH4,
    ele,
    ):
        """
        Computes the fluctuation of the reaction rate for the first reaction step in post-processing.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        prho : float
            Perturbation in density.
        pO2 : float
            Perturbation in O2 mass fraction.
        pCH4 : float
            Perturbation in CH4 mass fraction.
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Returns
        -------
        dQ1 : float
            Fluctuation of the reaction rate for reaction 1.
        """


        dT = -prho/MF['rho']*MF['T']
        return ((self.nu_O2_1+self.nu_CH4_1)*prho/MF['rho']\
                +self.beta1*dT/MF['T']\
                +self.Ta1*dT/MF['T']/MF['T']\
                +self.nu_O2_1*pO2/MF['O2']\
                +self.nu_CH4_1*pCH4/(MF['CH4']+self.epsilon))

    def postd_q2f_(self,
    MF,
    prho,
    pCO,
    pO2,
    ):
        """
        Computes the fluctuation of the forward reaction rate for the second reaction step in post-processing.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        prho : float
            Perturbation in density.
        pCO : float
            Perturbation in CO mass fraction.
        pO2 : float
            Perturbation in O2 mass fraction.

        Returns
        -------
        dQ2f : float
            Fluctuation of the forward reaction rate for reaction 2.
        """
        dT = -prho/MF['rho']*MF['T']
        return ((self.nu_O2_2+self.nu_CO_2)*prho/MF['rho']\
                +self.beta2*dT/MF['T']\
                +self.Ta2*dT/MF['T']/MF['T']\
                +self.nu_CO_2*pCO/(MF['CO'])\
                +self.nu_O2_2*pO2/MF['O2'])

    def postd_q2r_(self,
    MF,
    prho,
    pCO2,
    ):
        """
        Computes the fluctuation of the reverse reaction rate for the second reaction step in post-processing.

        Parameters
        ----------
        MF : dict
            Dictionary containing mean flow properties.
        prho : float
            Perturbation in density.
        pCO2 : float
            Perturbation in CO2 mass fraction.

        Returns
        -------
        dQ2r : float
            Fluctuation of the reverse reaction rate for reaction 2.
        """
        dK2rdT=self.d_k2rd_t_(self.T)
        dT = -prho/MF['rho']*MF['T']
        return (\
                +self.nu_CO2_2*pCO2/(MF['CO2'])\
                +(self.nu_CO2_2)*prho/MF['rho']\
                +dT*dK2rdT\
                )


    def test_m(self):
        """
        Prints the maximum reaction rates of the base flow for diagnostic purposes.

        Notes
        -----
        This method outputs the maximum values of reaction rates for reaction 1, reaction 2 forward, and reaction 2 inverse.
        """
        print('Reaction rates of base flow (maximum values)')
        print('reaction 1:')
        print(np.max(self.Q1.vector()[:]))
        print('reaction 2 forward:')
        print(np.max(self.Q2f.vector()[:]))
        print('reaction 2 inverse:')
        print(np.max(self.Q2r.vector()[:]))



#       #fitting functions
#       #BFER reaction 1
#       self.phi01=1.1
#       self.phi11=1.13
#       self.phi21=1.6
#       self.segma01=0.09
#       self.segma11=0.03
#       self.segma21=0.22
#       self.B1=0.37
#       self.C1=6.7

#       #BFER reaction 2
#       self.phi02=0.95
#       self.phi12=1.3
#       self.phi22=1.2
#       self.phi32=1.2
#       self.segma02=0.08
#       self.segma12=0.04
#       self.segma22=0.24
#       self.segma32=0.05
#       self.B2=2.5e-5
#       self.C2=8.7e-3

#to define the norm of heat release
# def addForm(self, MF, testf, trialf, solutionList):
#   self.i_rho=solutionList.index('rho')
#   self.i_CH4=solutionList.index('CH4')
#   self.i_O2=solutionList.index('O2')
#   self.i_CO2=solutionList.index('CO2')
#   self.i_CO=solutionList.index('CO')
#   dQ1=self.dQ1_(MF, trialf)
#   form = testf[self.i_rho]*dQ1*self.Q1*self.h1*dx

#   dQ2f=self.dQ2f_(MF, trialf)
#   dQ2r=self.dQ2r_(MF, trialf)
#   dQ2=dQ2f*self.Q2f-dQ2r*self.Q2r

#   form += testf[self.i_rho]*dQ2*self.h2*dx
#   return form
