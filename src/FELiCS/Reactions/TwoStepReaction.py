#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:31:52 2020

@author: cwang
Similar code structure as Max's c2sm2.py
"""

#TwoStep reaction methane-air combustion
#adapted for BFER or 2S_CH4_CM2 model
#CERFACS reference https://www.cerfacs.fr/cantera/mechanisms/meth.php


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
	def __init__(self, ModelName):
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

# 		#PEA fitting functions (homogeneously equal to 1 evaluated for lean premixed case.)
		#if required, formula and parameters given at the end of script)
# 		self.phi=None
# 		self.f1=None
# 		self.f2=None


		self.Y_CH4_limited=None

		self.epsilon=0.0002 #correction for base flow CH4


		self.reactionName='TwoStep'
		print('Initializing reaction '+self.reactionName+' '+ModelName)
		self.ReadReactionDict(ModelName)


	def ReadReactionDict(self, ModelName):
		fileMixture = open(self.mixtureDirectory,'r')
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

		fileSpecies = open(self.speciesDirectory,'r')
		speciesDictDict = eval(fileSpecies.read())
		fileSpecies.close()
		self.WCH4 = speciesDictDict['CH4']['mol_weight']
		self.WO2 = speciesDictDict['O2']['mol_weight']
		self.WCO2 = speciesDictDict['CO2']['mol_weight']
		self.WCO = speciesDictDict['CO']['mol_weight']
		self.WH2O = speciesDictDict['H2O']['mol_weight']
		self.WN2 = speciesDictDict['N2']['mol_weight']

	def computeMeanField(self,MF,ele):
		self.Y_CH4_limited = Expression("Y_CH4_ + epsilon_", Y_CH4_=MF['CH4'],epsilon_=self.epsilon, degree=2)

		ele.interpolate(Expression("A*exp(-Ta/T_) * pow(T_, beta)* pow(rho_, nu_CH4 + nu_O2) * pow(Y_CH4_ / W_CH4, nu_CH4) * pow(Y_O2_ / W_O2, nu_O2)",\
							A=self.A1, T_=MF['T'], Ta=self.Ta1, beta=self.beta1, rho_=MF['rho'], nu_CH4=self.nu_CH4_1, nu_O2=self.nu_O2_1, \
								Y_CH4_=self.Y_CH4_limited , W_CH4=self.WCH4, Y_O2_=MF['O2'], W_O2=self.WO2, degree=2))

		self.Q1 = ele
		self.K2f=project(Expression("A*exp(-Ta/T_) * pow(T_, beta)", A=self.A2, T_=MF['T'], Ta=self.Ta2, beta=self.beta2, degree=2), ele)
		self.Q2f=project(Expression("K2f_*pow(rho_, nu_CO + nu_O2) * pow(Y_CO_ / W_CO, nu_CO) * pow(Y_O2_ / W_O2, nu_O2)",\
							K2f_=self.K2f, rho_=MF['rho'], nu_CO=self.nu_CO_2, nu_O2=self.nu_O2_2, Y_CO_=MF['CO'],\
								W_CO=self.WCO, Y_O2_=MF['O2'], W_O2=self.WO2, degree=2), ele)

		#choose the base flow of temperature imported from AVBP. This can reproduce more precisely the base flow reaction rates
		self.T=MF['T']
		Tfield=MF['T'].vector[:]

		from Tabs.JANAFTabulate1971 import getEnthalpyPlusFormationJANAF1971,getEntropyJANAF1971


		lnexpInEqui_h = ((-np.array(list(getEnthalpyPlusFormationJANAF1971('CO2', Tfield)))/self.R/Tfield)*self.n_CO2_2\
				-(-np.array(list(getEnthalpyPlusFormationJANAF1971('CO', Tfield)))/self.R/Tfield)*self.n_CO_2\
				-(-np.array(list(getEnthalpyPlusFormationJANAF1971('O2', Tfield)))/self.R/Tfield)*self.n_O2_2)
		lnexpInEqui = (np.array(list(getEntropyJANAF1971('CO2',Tfield)))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971('CO2', Tfield)))/self.R/Tfield)*self.n_CO2_2\
				-(np.array(list(getEntropyJANAF1971('CO',Tfield)))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971('CO', Tfield)))/self.R/Tfield)*self.n_CO_2\
				-(np.array(list(getEntropyJANAF1971('O2',Tfield)))/self.R-np.array(list(getEnthalpyPlusFormationJANAF1971('O2', Tfield)))/self.R/Tfield)*self.n_O2_2

		nu_j=self.nu_CO2_2-self.nu_CO_2-self.nu_O2_2
		self.lnexpInEqui_h=Function(ele)
		self.lnexpInEqui_h.vector[:]=lnexpInEqui_h
		self.expInEqui=Function(ele)
		self.expInEqui.vector[:]=np.exp(lnexpInEqui)
		self.Equi=Function(ele)
		self.Equi.vector[:]=(self.p0/self.R/Tfield)**nu_j*np.exp(lnexpInEqui)

		self.K2r=project(Expression("k2f_/Equi_", k2f_=self.K2f, Equi_=self.Equi,  degree=2), ele)
		self.Q2r=project(Expression("K2r_*pow(rho_, nu_CO2) * pow(Y_CO2_ / W_CO2, nu_CO2) ",\
							K2r_=self.K2r, rho_=MF['rho'], nu_CO2=self.nu_CO2_2, Y_CO2_=MF['CO2'], W_CO2=self.WCO2, degree=2), ele)



		return self.Q1, self.Q2f, self.Q2r


	def addReaction(self, MF, testf, trialf, solutionList,ele):
		self.i_rho=solutionList.index('rho')
		self.i_CH4=solutionList.index('CH4')
		self.i_O2=solutionList.index('O2')
		self.i_CO2=solutionList.index('CO2')
		self.i_CO=solutionList.index('CO')

		#dQ1
		dQ1=self.dQ1_(MF, trialf)
		form = -testf[self.i_rho]*dQ1*self.Q1*self.h1*dx\
			-testf[self.i_CH4]*dQ1*self.Q1*self.n_CH4_1*self.WCH4*dx\
			-testf[self.i_O2]*dQ1*self.Q1*self.n_O2_1*self.WO2*dx\
			+testf[self.i_CO]*dQ1*self.Q1*self.n_CO_1*self.WCO*dx


# 		#dQ2f
		dQ2f=self.dQ2f_(MF, trialf)
		dQ2r=self.dQ2r_(MF, trialf)
		dQ2=dQ2f*self.Q2f-dQ2r*self.Q2r

		form += \
			-testf[self.i_O2]*dQ2*self.n_O2_2*self.WO2*dx\
			-testf[self.i_CO]*dQ2*self.n_CO_2*self.WCO*dx\
			+testf[self.i_CO2]*dQ2*self.n_CO2_2*self.WCO2*dx\
			-testf[self.i_rho]*dQ2*self.h2*dx
		return form

	def dQ1_(self, MF, trialf):
		dT = -trialf['rho']/MF['rho']*MF['T']
		return ((self.nu_O2_1+self.nu_CH4_1)*trialf['rho']/MF['rho']\
				+self.beta1*dT/MF['T']\
				+self.Ta1*dT/MF['T']/MF['T']\
				+self.nu_O2_1*trialf['O2']/MF['O2']\
				+self.nu_CH4_1*trialf['CH4']/(MF['CH4']+self.epsilon))#((MF['O2']-0.0447)/1.78e-1*4.45e-2))#/MF['CH4'])#self.Y_CH4_limited)


	def dQ2f_(self, MF, trialf):
		dT = -trialf['rho']/MF['rho']*MF['T']
		return ((self.nu_O2_2+self.nu_CO_2)*trialf['rho']/MF['rho']\
				+self.beta2*dT/MF['T']\
				+self.Ta2*dT/MF['T']/MF['T']\
				+self.nu_CO_2*trialf['CO']/(MF['CO'])\
				+self.nu_O2_2*trialf['O2']/MF['O2'])

	def dQ2r_(self, MF, trialf):
		dK2rdT=self.dK2rdT_(self.T)
		dT = -trialf['rho']/MF['rho']*MF['T']
		return (\
				+self.nu_CO2_2*trialf['CO2']/(MF['CO2'])\
				+self.nu_CO2_2*trialf['rho']/MF['rho']\
				+dT*dK2rdT\
				)
# 		#if we neglect dK2rdT
# 		return (\
# 				+self.nu_CO2_2*trialf['CO2']/(MF['CO2'])\
# 				+self.nu_CO2_2*trialf['rho']/MF['rho']\
# 				+dT*dK2rdT\
# 				)

	def dK2fdT_(self,T):
		return Expression("(Ta2 / pow(T_, 2)+beta2/T_)", Ta2=self.Ta2, T_=T, beta2=self.beta2, degree=2)

	def dK2rdT_(self, T):

		return Expression("Kf2dT_ -KcdT_/Kc",Kf2dT_=self.dK2fdT_(T), Kc=self.Equi, K2f_=self.K2f, K2r_=self.K2r,\
 					KcdT_=self.dEquidT_(self.T),  degree=2)
#		#if we neglect KcdT
#		return Expression("Kf2dT_ ",Kf2dT_=self.dK2fdT_(T), Kc=self.Equi, K2f_=self.K2f, K2r_=self.K2r,\
#					KcdT_=self.dEquidT_(self.T),  degree=2)

	def dEquidT_(self, T):
 		return Expression(
 		"(Kp2dT - Kp2/T_ * (n_CO2_2 - n_O2_2 - n_CO_2)) * pow(pa / (R * T_), n_CO2_2-n_O2_2-n_CO_2)",
 		Kp2dT=self.dexpInEquidT_(T), Kp2=self.expInEqui, pa=self.p0, R=self.R, T_=T,
 		n_O2_2=self.n_O2_2, n_CO_2=self.n_CO_2, n_CO2_2=self.n_CO2_2,
 		degree=2)
#		#if we neglect Kp2dT
#		return Expression(
#		"(- Kp2/T_ * (n_CO2_2 - n_O2_2 - n_CO_2)) * pow(pa / (R * T_), n_CO2_2-n_O2_2-n_CO_2)",
#		 Kp2=self.expInEqui, pa=self.p0, R=self.R, T_=T,
#		n_O2_2=self.n_O2_2, n_CO_2=self.n_CO_2, n_CO2_2=self.n_CO2_2,
#		degree=2)

	def dexpInEquidT_(self, T):
		return Expression("- Kp2 /  T_ * H_",
		Kp2=self.expInEqui, R=self.R, T_=T,
		H_=self.lnexpInEqui_h,
		degree=2)


	def postHeatRelaese(self, MF, prho, pCH4, pO2, pCO, pCO2, ele):
		#evalute fluctuation HeatRelease in post-processing
		#dQ1
		dQ1=self.postdQ1_(MF, prho, pO2, pCH4, ele)
		form = dQ1*self.Q1*self.h1

		#dQ2
		dQ2f=self.postdQ2f_(MF, prho, pCO, pO2)
		dQ2r=self.postdQ2r_(MF, prho, pCO2)
		form += (dQ2f*self.Q2f-dQ2r*self.Q2r)*self.h2

		return project(form,ele)

	def postdQ1_(self, MF, prho, pO2, pCH4, ele):


		dT = -prho/MF['rho']*MF['T']
		return ((self.nu_O2_1+self.nu_CH4_1)*prho/MF['rho']\
				+self.beta1*dT/MF['T']\
				+self.Ta1*dT/MF['T']/MF['T']\
				+self.nu_O2_1*pO2/MF['O2']\
				+self.nu_CH4_1*pCH4/(MF['CH4']+self.epsilon))

	def postdQ2f_(self, MF, prho, pCO, pO2):
		dT = -prho/MF['rho']*MF['T']
		return ((self.nu_O2_2+self.nu_CO_2)*prho/MF['rho']\
				+self.beta2*dT/MF['T']\
				+self.Ta2*dT/MF['T']/MF['T']\
				+self.nu_CO_2*pCO/(MF['CO'])\
				+self.nu_O2_2*pO2/MF['O2'])

	def postdQ2r_(self, MF, prho, pCO2):
		dK2rdT=self.dK2rdT_(self.T)
		dT = -prho/MF['rho']*MF['T']
		return (\
				+self.nu_CO2_2*pCO2/(MF['CO2'])\
				+(self.nu_CO2_2)*prho/MF['rho']\
				+dT*dK2rdT\
				)


	def testM(self):
		print('Reaction rates of base flow (maximum values)')
		print('reaction 1:')
		print(np.max(self.Q1.vector()[:]))
		print('reaction 2 forward:')
		print(np.max(self.Q2f.vector()[:]))
		print('reaction 2 inverse:')
		print(np.max(self.Q2r.vector()[:]))



# 		#fitting functions
# 		#BFER reaction 1
# 		self.phi01=1.1
# 		self.phi11=1.13
# 		self.phi21=1.6
# 		self.segma01=0.09
# 		self.segma11=0.03
# 		self.segma21=0.22
# 		self.B1=0.37
# 		self.C1=6.7

# 		#BFER reaction 2
# 		self.phi02=0.95
# 		self.phi12=1.3
# 		self.phi22=1.2
# 		self.phi32=1.2
# 		self.segma02=0.08
# 		self.segma12=0.04
# 		self.segma22=0.24
# 		self.segma32=0.05
# 		self.B2=2.5e-5
# 		self.C2=8.7e-3

#to define the norm of heat release
# def addForm(self, MF, testf, trialf, solutionList):
# 	self.i_rho=solutionList.index('rho')
# 	self.i_CH4=solutionList.index('CH4')
# 	self.i_O2=solutionList.index('O2')
# 	self.i_CO2=solutionList.index('CO2')
# 	self.i_CO=solutionList.index('CO')
# 	dQ1=self.dQ1_(MF, trialf)
# 	form = testf[self.i_rho]*dQ1*self.Q1*self.h1*dx

# 	dQ2f=self.dQ2f_(MF, trialf)
# 	dQ2r=self.dQ2r_(MF, trialf)
# 	dQ2=dQ2f*self.Q2f-dQ2r*self.Q2r

# 	form += testf[self.i_rho]*dQ2*self.h2*dx
# 	return form
