#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from dolfinx.fem import (
    Function,
    Expression,
)
from ufl import (
    dx,
    exp,
    conj,
)


class GlobalReaction():
    """
    Represents a global one-step chemical reaction mechanism.

    This class encapsulates the parameters and methods required to model a global reaction, including reading reaction data, computing mean fields, and handling heat release.

    **Initialize the GlobalReaction object**

    Parameters
    ----------
    reaction_mechanism : dict
        Dictionary containing reaction mechanism parameters.

    Attributes
    ----------
    mixtureDirectory : str
        Directory name for mixture data.
    speciesDirectory : str
        Directory name for species data.
    A : float or None
        Pre-exponential factor for the reaction.
    Ta : float or None
        Activation temperature.
    a : float or None
        Stoichiometric coefficient for oxygen.
    b : float or None
        Stoichiometric coefficient for carbon.
    beta : float or None
        Temperature exponent.
    h0 : float or None
        Heat of reaction.
    st_C : float or None
        Stoichiometric coefficient for carbon species.
    st_O : float or None
        Stoichiometric coefficient for oxygen species.
    WO : float or None
        Molecular weight of oxygen.
    WC : float or None
        Molecular weight of carbon.
    reactionName : str
        Name of the reaction.
    Q : dolfinx.fem.Function or None
        Mean field reaction rate.
    i_rho : int or None
        Index of density in solution list.
    i_C : int or None
        Index of CH4 in solution list.
    """

    def __init__(self, reaction_mechanism):
        """
        Initializes the GlobalReaction instance.

        Parameters
        ----------
        reaction_mechanism : dict
            Dictionary containing reaction mechanism parameters.
        """

        self.mixtureDirectory = "Mixture" #to be put in param
        self.speciesDirectory = "Species" #to be put in param

        self.A=None
        self.Ta=None
        self.a=None
        self.b=None
        self.beta=None
        self.h0=None
        self.st_C=None
        self.st_O=None
        self.WO=None
        self.WC=None

        self.reactionName='Global Reaction'
        print('Initializing reaction '+self.reactionName)
        self.ReadReactionDict(reaction_mechanism)

        self.Q=None
        self.i_rho=None
        self.i_C=None

    def ReadReactionDict(self, reaction):
        """
        Reads reaction parameters from a dictionary and updates class attributes.

        Parameters
        ----------
        reaction : dict
            Dictionary containing reaction parameters.

        Notes
        -----
        This method sets the reaction coefficients and molecular weights based on the provided reaction dictionary.
        """

        self.A=reaction['reac_preexp']
        self.Ta=reaction['reac_act_tem']
        self.a=reaction['reac_nu_O']
        self.b=reaction['reac_nu_C']
        self.beta=reaction['reac_exp_tem']
        self.st_O=reaction['reac_st_O']
        self.st_C=reaction['reac_st_C']
        self.h0=reaction['reac_h0']
        educt_C=reaction['educt_C']
        educt_O=reaction['educt_O']

        fileSpecies = open(self.speciesDirectory,'r')
        speciesDictDict = eval(fileSpecies.read())
        fileSpecies.close()
        self.WO=speciesDictDict[educt_O]['mol_weight']
        self.WC=speciesDictDict[educt_C]['mol_weight']

    def computeMeanField(self,mean,ele):
        """
        Computes the mean field reaction rate and interpolates it as a function.

        Parameters
        ----------
        mean : object
            Object containing mean flow properties (e.g., temperature, density, species mass fractions).
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Returns
        -------
        Q : dolfinx.fem.Function
            Interpolated mean field reaction rate.

        Notes
        -----
        Uses the reaction parameters and mean flow properties to construct the reaction rate expression.
        """

        self.Q = Function(ele)
        # https://jorgensd.github.io/dolfinx-tutorial/chapter1/membrane_code.html#interpolation-of-a-ufl-expression
        expressionUFL = self.A*exp(-self.Ta/mean.T) * mean.rho**(self.b + self.a) * (mean.Y('CH4') / self.WC) **self.b * (mean.Y('O2') / self.WO)**self.a
        expr = Expression(expressionUFL, ele.element.interpolation_points())
        self.Q.interpolate(expr)
        #self.Q = ele
        return self.Q

    def addReaction(self, mean, testf, fluc, solutionList):
        """
        Adds the reaction term to the weak form.

        Parameters
        ----------
        mean : object
            Object containing mean flow properties.
        testf : list
            List of test functions.
        fluc : object
            Object containing fluctuation properties.
        solutionList : list of str
            List of solution variable names.

        Returns
        -------
        reaction_term : ufl.Form
            Weak form expression for the reaction term.

        Notes
        -----
        The method computes the indices for relevant variables and constructs the weak form using the reaction rate and heat release.
        """

        self.i_rho=solutionList.index('rho')
        self.i_C=solutionList.index('CH4')
        dQ=self.dQ_(mean, fluc)
        return -conj(testf[self.i_rho])*(dQ*self.Q)*self.h0*dx-conj(testf[self.i_C])*(dQ*self.Q)*self.st_C*self.WC*dx

    def dQ_(self, mean, fluc):
        """
        Computes the fluctuation of the reaction rate.

        Parameters
        ----------
        mean : object
            Object containing mean flow properties.
        fluc : object
            Object containing fluctuation properties.

        Returns
        -------
        dQ : float
            Fluctuation of the reaction rate.

        Notes
        -----
        The calculation uses stoichiometric coefficients and mean/fluctuation values for temperature and species.
        """

        dO2_= fluc.Y('CH4')/(self.st_C*self.WC)*self.st_O*self.WO
        dT_ = -fluc.rho/mean.rho*mean.T
        return ((self.a+self.b)*fluc.rho/mean.rho\
                +self.beta*dT_/mean.T\
                +self.Ta*dT_/mean.T/mean.T\
                +self.a*dO2_/mean.Y('O2')\
                +self.b*fluc.Y('CH4')/mean.Y('CH4'))

    def postHeatRelease(self, mean, prho, pCH4, ele):
        """
        Computes the post-processed heat release form.

        Parameters
        ----------
        mean : object
            Object containing mean flow properties.
        prho : float
            Perturbation in density.
        pCH4 : float
            Perturbation in CH4 mass fraction.
        ele : dolfinx.fem.FunctionSpace
            Function space for interpolation.

        Notes
        -----
        The `ele` parameter is included for compatibility with the original method signature; however, it is not used in this mehtod.

        Returns
        -------
        heat_release_field : dolfinx.fem.Function
            Interpolated post-processed heat release form.
        """

        dQ=self.postdQ(mean, prho, pCH4)
        form =-(dQ*self.Q)*self.h0
        return dQ.interpolate(form)

    def postdQ(self, mean, prho, pCH4):
        """
        Computes the post-processed fluctuation of the reaction rate.

        Parameters
        ----------
        mean : object
            Object containing mean flow properties.
        prho : float
            Perturbation in density.
        pCH4 : float
            Perturbation in CH4 mass fraction.

        Returns
        -------
        dQ : float
            Post-processed fluctuation of the reaction rate.
        """

        dO2_= pCH4/(self.st_C*self.WC)*self.st_O*self.WO
        dT_ = -prho/mean.rho*mean.T
        return ((self.a+self.b)*prho/mean.rho\
                +self.beta*dT_/mean.T\
                +self.Ta*dT_/mean.T/mean.T\
                +self.a*dO2_/mean.Y('O2')\
                +self.b*pCH4/mean.Y('CH4'))

#CERFACS 1S_CH4_MP1 https://www.cerfacs.fr/cantera/mechanisms/meth.php
#       self.A=1.1e7
#       self.Ta=20000/1.987
#       self.a=0.5
#       self.b=1.0
#       self.A=6.7e9
#       self.Ta=48400/1.987
#       self.a=1.3
#       self.b=0.2
#       self.beta=0.0
#       self.h0=-2*238.922e3-393.151e3+66.911e3#+74.8e3
#       self.st_C=1.0
#       self.st_O=2.0
