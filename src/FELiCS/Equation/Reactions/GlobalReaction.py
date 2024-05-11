#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 15:39:58 2020

@author: cwang
"""

#OneStep global reaction
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
    def __init__(self, reaction_mechanism):
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
        #fileMixture = open('Mixture','r')
        #mixtureDictDict = eval(fileMixture.read())
        #fileMixture.close()
        #ModelName='WestbrookDryer_Max'
        #MD = mixtureDictDict[ModelName]
        #print(MD)
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
        #import pdb
        #import ufl
        self.Q = Function(ele)
        # https://jorgensd.github.io/dolfinx-tutorial/chapter1/membrane_code.html#interpolation-of-a-ufl-expression
        expressionUFL = self.A*exp(-self.Ta/mean.T) * mean.rho**(self.b + self.a) * (mean.Y('CH4') / self.WC) **self.b * (mean.Y('O2') / self.WO)**self.a
        expr = Expression(expressionUFL, ele.element.interpolation_points())
        self.Q.interpolate(expr)
        #self.Q = ele
        return self.Q

    def addReaction(self, mean, testf, fluc, solutionList):
        self.i_rho=solutionList.index('rho')
        self.i_C=solutionList.index('CH4')
        dQ=self.dQ_(mean, fluc)
        return -conj(testf[self.i_rho])*(dQ*self.Q)*self.h0*dx-conj(testf[self.i_C])*(dQ*self.Q)*self.st_C*self.WC*dx

    def dQ_(self, mean, fluc):
        dO2_= fluc.Y('CH4')/(self.st_C*self.WC)*self.st_O*self.WO
        dT_ = -fluc.rho/mean.rho*mean.T
        return ((self.a+self.b)*fluc.rho/mean.rho\
                +self.beta*dT_/mean.T\
                +self.Ta*dT_/mean.T/mean.T\
                +self.a*dO2_/mean.Y('O2')\
                +self.b*fluc.Y('CH4')/mean.Y('CH4'))

    def postHeatRelease(self, mean, prho, pCH4, ele):
        dQ=self.postdQ(mean, prho, pCH4)
        form =-(dQ*self.Q)*self.h0
        return dQ.interpolate(form)

    def postdQ(self, mean, prho, pCH4):
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
