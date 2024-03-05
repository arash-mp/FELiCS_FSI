#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:40:33 2021

@author: cwang
"""

import numpy as np
#JANAF polynomial (Max's flame OpenFOAM)
def getJANAFPolynomialTable():
	#Ref temperature T=0 K
	#cp=R*(a1+a2*T+a3*T**2+a4*T**3+a5*T**4)
	#h=R*(a1*T+a2*T**2/2+a3*T**3/3+a4*T**4/4+a5*T**5/5+a6)
	#s=R*(a1*ln(T)+a2*T+a3*T**2/2+a4*T**3/3+a5*T**4/4+a7)
	#for T>1000K [a1, a2, a3, a4, a5, a6, a7\
	#for T<1000K a1, a2, a3, a4, a5, a6, a7]
	tempDict={}
	tempDict['O2']=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00,\
		 3.78245636E+00,-2.99673416E-03, 9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00,1000, 0.0319988e0]
	tempDict['H2O']=[3.03399249E+00, 2.17691804E-03,-1.64072518E-07,-9.70419870E-11, 1.68200992E-14,-3.00042971E+04,4.96677010E+00,\
		 4.19864056E+00,-2.03643410E-03, 6.52040211E-06,-5.48797062E-09, 1.77197817E-12,-3.02937267E+04,-8.49032208E-01,1000, 3.002628e-02]

	tempDict['CH4']=[7.48514950E-02, 1.33909467E-02,-5.73285809E-06, 1.22292535E-09,-1.01815230E-13,-9.46834459E+03,1.84373180E+01,\
				  5.14987613E+00,-1.36709788E-02, 4.91800599E-05,-4.84743026E-08, 1.66693956E-11,-1.02466476E+04,-4.64130376E+00,1000, 0.0160423e0]
	tempDict['CO']=[2.71518561E+00, 2.06252743E-03,-9.98825771E-07, 2.30053008E-10,-2.03647716E-14,-1.41518724E+04,7.81868772E+00,\
		 3.57953347E+00,-6.10353680E-04, 1.01681433E-06,9.07005884E-10,-9.04424499E-13,-1.43440860E+04,3.50840928E+00,1000, 0.0280140e0]
	tempDict['CO2']=[3.85746029E+00, 4.41437026E-03,-2.21481404E-06, 5.23490188E-10,-4.72084164E-14,-4.87591660E+04, 2.27163806E+00,\
		 2.35677352E+00, 8.98459677E-03,-7.12356269E-06,2.45919022E-09,-1.43699548E-13,-4.83719697E+04,9.90105222E+00,1000, 0.0440098e0]
	tempDict['N2']=[0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.06753351E-13, -0.09227977E+04, 0.05980528E+02,\
		0.03298677E+02, 0.14082404E-02, -0.03963222E-04, 0.05641515E-07, -0.02444854E-10, -0.10208999E+04, 0.03950372E+02,1000, 0.0280134e0]
	return tempDict

def getEnthalpyJANAFelement(coeffs, T):
	R=8.31446261815324
	if T>coeffs[-2]:
		[a1, a2, a3, a4, a5, a6]=coeffs[0:6]
	else:
		[a1, a2, a3, a4, a5, a6]=coeffs[7:13]
	return (R*(a1*T+a2*T**2/2+a3*T**3/3+a4*T**4/4+a5*T**5/5+a6)-R*coeffs[12])

def getEnthalpyJANAF(species, T):
	coeffs=getJANAFPolynomialTable()[species]
	return map(lambda x: getEnthalpyJANAFelement(coeffs,x),T)



def getCpJANAFelement(coeffs, T):
	R=8.31446261815324
	if T>coeffs[-2]:
		[a1, a2, a3, a4, a5, a6]=coeffs[0:6]
	else:
		[a1, a2, a3, a4, a5, a6]=coeffs[7:13]
	return R*(a1+a2*T+a3*T**2+a4*T**3+a5*T**4)/coeffs[-1]

def getCpJANAF(species, T):
	coeffs=getJANAFPolynomialTable()[species]
	return map(lambda x: getCpJANAFelement(coeffs,x),T)

def getEntropyJANAFelement(coeffs, T):
	R=8.31446261815324
	if T>coeffs[-2]:
		[a1, a2, a3, a4, a5]=coeffs[0:5]
		a7=coeffs[6]
	else:
		[a1, a2, a3, a4, a5]=coeffs[7:12]
		a7=coeffs[13]
	return R*(a1*np.log(T)+a2*T+a3*T**2/2+a4*T**3/3+a5*T**4/4+a7)

def getEntropyJANAF(species, T):
	coeffs=getJANAFPolynomialTable()[species]
	return np.array(list(map(lambda x: getEntropyJANAFelement(coeffs,x),T)))

#test
# import numpy as np
# #T=np.array([300,400])
# T=np.linspace(0,5000,51)
# T[-1]=4999
# h=np.array(list(getEnthalpyJANAF('CH4', T)))#-np.array(list(getEnthalpyJANAF('CH4', T0)))
# cp=np.array(list(getCpJANAF('CH4', T)))
# ss=np.array(list(getEntropyJANAF('CH4', T)))
# cp1971=np.array(list(getCpJANAF1971('CH4', T)))
# h1971=np.array(list(getEnthalpyJANAF1971('CH4', T)))
# ss1971=np.array(list(getEntropyJANAF1971('CH4', T)))

# hCH4=[0.0000000000e+00,  0.3326000000e+04,  0.6656000000e+04,  0.1009000000e+05,  0.1388500000e+05,  0.1822400000e+05,  0.2315400000e+05,  0.2865900000e+05,  0.3469900000e+05,  0.4122900000e+05,  0.4820300000e+05,  0.5557300000e+05,  0.6329400000e+05,  0.7132600000e+05,  0.7963200000e+05,  0.8817700000e+05,  0.9693400000e+05,  0.1058770000e+06,  0.1149840000e+06,  0.1242360000e+06,  0.1336160000e+06,  0.1431110000e+06,  0.1527080000e+06,  0.1623950000e+06,  0.1721650000e+06,  0.1820080000e+06,  0.1919170000e+06,  0.2018860000e+06,  0.2119090000e+06,  0.2219820000e+06,  0.2321000000e+06,  0.2422590000e+06,  0.2524550000e+06,  0.2626860000e+06,  0.2729490000e+06,  0.2832410000e+06,  0.2935600000e+06,  0.3039050000e+06,  0.3142720000e+06,  0.3246610000e+06,  0.3350690000e+06,  0.3454970000e+06,  0.3559420000e+06,  0.3664030000e+06,  0.3768790000e+06,  0.3873690000e+06,  0.3978730000e+06,  0.4083890000e+06,  0.4189170000e+06,  0.4294560000e+06,  0.4400060000e+06]
# sCH4=[000.000e0,149.500e0,172.577e0,186.472e0,197.356e0, 207.014e0,215.987e0,224.461e0,232.518e0,240.205e0,\
# 				 247.549e0,254.570e0,261.287e0,267.714e0,273.868e0, 279.763e0,285.413e0,290.834e0,296.039e0,301.041e0,\
# 					305.853e0,310.485e0,314.949e0,319.255e0,323.413e0, 327.431e0,331.317e0,335.080e0,338.725e0,342.260e0,\
# 					345.690e0,349.021e0,352.258e0,355.406e0,358.470e0, 361.453e0,364.360e0,367.194e0,369.959e0,372.658e0,\
# 					375.293e0,377.868e0,380.385e0,382.846e0,385.255e0, 387.612e0,389.921e0,392.182e0,394.399e0,396.572e0, 398.703e0]



# import matplotlib.pyplot as plt
# plt.figure()
# plt.plot(T,h, label='polynomial')
# plt.plot(T,h1971, label='1971')
# plt.plot(T,hCH4, label='Thomas table')
# plt.legend()
# plt.xlabel('temperature')
# plt.ylabel('h')
# plt.savefig('JANAF_enthalpy')

# plt.figure()
# plt.plot(T,cp, label='polynomial')
# plt.plot(T,cp1971, label='1971')
# plt.legend()
# plt.xlabel('temperature')
# plt.ylabel('cp')
# plt.savefig('JANAF_Cp')

# plt.figure()
# plt.plot(T,ss, label='polynomial')
# plt.plot(T,sCH4, label='Thomas table')
# plt.plot(T,ss1971, label='1971')
# plt.legend()
# plt.xlabel('temperature')
# plt.ylabel('s')
# plt.savefig('JANAF_entropy')