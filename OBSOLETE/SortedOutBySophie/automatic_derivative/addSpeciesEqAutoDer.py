#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 15:45:55 2020

@author: cwang
"""

from fenics import as_vector,inner,derivative, grad, dx
from ufl import i

def addSpeciesEq(self,fluc,X,mean,param, rup):
	#nonlinear space: rup
	ru  = as_vector((rup[0], rup[1]))
	rrho = rup[3]
	
	#derived variables
	print('Sutherland law for species perturbation')
	Ts=170.672
	perturbT=(mean['T']+3*Ts)/(2*(mean['T']+Ts))*(-fluc['rho']/mean['rho']) #visco fluctuation factor(Sutherland)
	
	for specie_name in param.SpeciesList:
		i_eqn=param.SolutionList.index(specie_name)
		ry=rup[i_eqn+1] #trial function
		y=X[i_eqn] #testfunction

		#convection
		#F = X[i_eqn]*rrho*inner(ru,grad(ry))*dx 
		#F = -y*rrho*inner(ru,grad(fluc[specie_name]))*dx 
		#diffusion
		Sc=0.7
		#F = -mean['MuTot']/Sc*(inner(grad(ry),grad(y)))*dx
		#self.A_imag_vf += derivative(F, rup)
		
		#time derivative
		##self.B_real_vf += X[i_eqn]*fluc[specie_name]*self.rho*dx
		#viscosity fluctuation
		#self.A_imag_vf += -mean['MuTot']/Sc*perturbT*(inner(grad(mean[specie_name]),grad(X[i_eqn])))*dx

