#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 12:30:35 2019

@author: cwang, T. Kaiser
"""
import Mesh
import Import
import DefineFEMSpaces
from ExportSolution import *
from param_BCs import *
import global_variables as glob
import parameters as param
import functions


from dolfin import *
import pyvtk
import numpy as np
#import pickle
from numpy.random import seed, random
import matplotlib.pyplot as plt




print("LoadingMesh...")
mesh=Mesh("BunsenFlame.xml")






# Define Boundary conditions
def symmetry(x, on_boundary):
	return x[1] < DOLFIN_EPS and on_boundary
def wall(x, on_boundary):
	return x[0] < DOLFIN_EPS and x[1]> 0.006-DOLFIN_EPS and on_boundary
def inlet(x, on_boundary):
	return x[0] < (-0.01+DOLFIN_EPS) and on_boundary
def outlet(x, on_boundary):
	return x[0] > (0.09-DOLFIN_EPS) and on_boundary
def up_bc(x, on_boundary):
	return  x[1]>0.0125-DOLFIN_EPS and on_boundary

#print Meanflow
def plotMF(MeanFlowDict):
	nMeanFlow=len(MeanFlowDict.keys())
	#print(nMeanFlow)
	fig = plt.figure()
	plt.set_cmap('coolwarm')
	if (param.PlotFlows):
		i=0
		for name in MeanFlowDict.keys():
			#name=nameListMean[i]
			#name=list(MeanFlowDict.keys())[i]
			ax = fig.add_subplot(nMeanFlow,1,i+1)
			i=i+1
	#		print(name)
			cs=plot(MeanFlowDict[name])


			cbar = fig.colorbar(cs,ticks=[np.min(MeanFlowDict[name].vector()[:]), np.max(MeanFlowDict[name].vector()[:])])
			plt.title(name)
		plt.show()


def temperaTure():
	p_constant=101300
	return p_constant/rho/(8314.4598/28.949)

#viscosity
def visco():
	C1=1.672e-6
	C2=170.7
	return C1*temperaTure()**(3/2)/(temperaTure()+C2)*relax_visco

#increment of visco
def viscop():
	C2=170.7
	return -rhop/rho*visco()*(temperaTure()+3*C2)/(2*temperaTure()+2*C2)

#reaction rate with OneStep model
def omegaReaction():
	ReactionField=Function(FEMSpaces.P1)
	ReactionField=Reaction['A']\
		*rho**(Reaction['a']+Reaction['b'])\
		*YO2**Reaction['a']\
		*YCH4**Reaction['b']\
		/MolecularMass[Reaction['educts'][0]]**Reaction['a']\
		/MolecularMass[Reaction['educts'][1]]**Reaction['b']\
		*exp(-Reaction['E_a']/(1.013e5/rho/(8314.4598/28.949))/Reaction['R'])
	return ReactionField




	# Define function spaces (P2-P1)
V = VectorElement("CG", mesh.ufl_cell(), 2)
R= FiniteElement("CG", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([V,R,R,R,R]))


FEMSpaces=DefineFEMSpaces.FEMSpacesClass(param,mesh)
#plotMF(MF)




base=Function(W)
VF=FunctionSpace(mesh,V)
RF=FunctionSpace(mesh,R)

#Two options : Assign start solution OR Read solution

###Assign start solution
#
##Assign velocity
#u_init=Function(VF)
##Assign seperately ux and uy profile
##ux profile
#UUFS=FunctionSpace(mesh,FiniteElement("CG", mesh.ufl_cell(), 2))
#UUFx=Function(UUFS)
##ppp.interpolate(Expression('(x[1]<0.006+DOLFIN_EPS) ? 2*0.5*(tanh(2.5*(1*0.006/x[1]-x[1]/0.006))) : 0.0 ',degree=2)) #jet profile
#UUFx.interpolate(Expression('(x[1]<0.006+DOLFIN_EPS) ? (1-x[1]/0.006*x[1]/0.006) : 0.0 ',degree=2)) #parabolic profile
##uy profile
#UUFy=Function(UUFS)
#UUFyy.interpolate(Expression('0.0',degree=2))
#assign(u_init,[ppp,pppy])
##Assign velocity field with AVBP data
##assign(u_init,[MF['ux'],pppy])
#
#
##Assign pressure
#p_init=Function(RF)
##assign(p_init,MF['p'])
#p_init.interpolate(Expression('0.0',degree=1))
#
#
##Assign density
#rho_init=Function(RF)
##rho_init.interpolate(Expression('(x[0]<-DOLFIN_EPS)? 1.1530 : 0.1676',degree=1))
##rho_init.interpolate(Expression('(x[1]<0.006+DOLFIN_EPS) ? 1.153 : 0.8',degree=1))
#assign(rho_init,MF['rho'])
#
##Assign YO2
#YO2_init=Function(RF)
##YO2_init.interpolate(Expression('(x[0]<-DOLFIN_EPS) ? 0.2232 : 0.0000001',degree=1))
#assign(YO2_init,MF['O2'])
#
##Assign CH4
#YCH4_init=Function(RF)
##YCH4_init.interpolate(Expression('(x[0]<-DOLFIN_EPS)  ? 0.04472 : 0.0000001',degree=1))
#assign(YCH4_init,MF['CH4'])


#Read solution


str_E_Re='Sun1_h1_E0p76_i1'
u_init=Function(VF,str_E_Re+'_u.xml')
p_init=Function(RF,str_E_Re+'_p.xml')
rho_init=Function(RF,str_E_Re+'_rho.xml')
YCH4_init=Function(RF,str_E_Re+'_yCH4.xml')
YO2_init=Function(RF,str_E_Re+'_yO2.xml')

assigner = FunctionAssigner(W, [VF,RF,RF,RF,RF])
assigner.assign(base,[u_init,p_init,rho_init,YO2_init,YCH4_init])
u=as_vector((base[0],base[1]))
p=base[2]
rho=base[3]
YO2=base[4]
YCH4=base[5]


# Define trial and test functions
perturbation = TrialFunction(W)
up=as_vector((perturbation[0],perturbation[1]))
pp=perturbation[2]
rhop=perturbation[3]
YO2p=perturbation[4]
YCH4p=perturbation[5]

testMixed = TestFunction(W)
v=as_vector((testMixed[0],testMixed[1]))
q=testMixed[2]
s=testMixed[3]
y1=testMixed[4]
y2=testMixed[5]


#plot init flow

fig, ax = plt.subplots()
plt.set_cmap('coolwarm')
plt.subplot(611)
cs=plot(base[0])
cbar = fig.colorbar(cs)
plt.subplot(612)
cs=plot(base[1])
cbar = fig.colorbar(cs)
plt.subplot(613)
cs=plot(base[2])
cbar = fig.colorbar(cs)
plt.subplot(614)
cs=plot(base[3])
cbar = fig.colorbar(cs)
plt.subplot(615)
cs=plot(base[4])
cbar = fig.colorbar(cs)
plt.subplot(616)
cs=plot(base[5])
cbar = fig.colorbar(cs)

# Define boundary conditions
bcs=[]
bcs.append(DirichletBC(W.sub(0), (0.0,0.0), wall ))
bcs.append(DirichletBC(W.sub(0), (0.0,0.0), inlet ))
bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, up_bc ))
bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, symmetry))

bcs.append(DirichletBC(W.sub(1), 0.0, inlet ))
bcs.append(DirichletBC(W.sub(1), 0.0, wall ))
bcs.append(DirichletBC(W.sub(1), 0.0, up_bc ))

bcs.append(DirichletBC(W.sub(2), 0.0, inlet ))
bcs.append(DirichletBC(W.sub(2), 0.0, wall ))
bcs.append(DirichletBC(W.sub(2), 0.0, up_bc ))


bcs.append(DirichletBC(W.sub(3), 0.0, inlet ))
bcs.append(DirichletBC(W.sub(4), 0.0, inlet ))
bcs.append(DirichletBC(W.sub(3), 0.0, wall ))
bcs.append(DirichletBC(W.sub(4), 0.0, wall ))
bcs.append(DirichletBC(W.sub(3), 0.0, up_bc ))
bcs.append(DirichletBC(W.sub(4), 0.0, up_bc ))




#parameters


tol=1e-2


Sc=0.7
Pr=0.7

relax_visco_obj=1  #See function visco()
relax_visco=relax_visco_obj

Ea_obj=48400
Reaction['E_a']=Ea_obj*0.755

A_obj=(6.7e9)
Reaction['A']=A_obj

kappa=1.4
h0=(817e3)
patom=1.013e5

io=0

while io<1:

	i=0
	err=1
	io=io+1

	while (err>tol and i<20):




	#Full reaction equation

	#Complete equation : NewtonNS
		NewtonNS=rho*inner(grad(u)*u,v)*dx+visco()*inner(grad(u),grad(v))*dx+inner(grad(p),v)*dx+1/3*visco()*(v[0].dx(0)*u[0].dx(0)+v[0].dx(0)*u[1].dx(1)+v[1].dx(1)*u[0].dx(0)+v[1].dx(1)*u[1].dx(1))*dx

		NewtonNS=NewtonNS-q*rho**2*inner(u,grad(rho))*dx
		NewtonNS=NewtonNS-1/Pr*visco()*(q*2*inner(grad(rho),grad(rho))+rho*inner(grad(rho),grad(q)))*dx
		NewtonNS=NewtonNS-q*rho**3*(kappa-1)/kappa/patom*h0*omegaReaction()*dx

		NewtonNS=NewtonNS+s*rho**3*div(u)*dx
		NewtonNS=NewtonNS-1/Pr*visco()*(s*2*inner(grad(rho),grad(rho))+rho*inner(grad(rho),grad(s)))*dx
		NewtonNS=NewtonNS-s*rho**3*(kappa-1)/kappa/patom*h0*omegaReaction()*dx

		NewtonNS=NewtonNS+y1*rho*inner(u,grad(YO2))*dx+visco()/Sc*inner(grad(YO2),grad(y1))*dx+y1*omegaReaction()*MolecularMass['O2']*2*dx
		NewtonNS=NewtonNS+y2*rho*inner(u,grad(YCH4))*dx+visco()/Sc*inner(grad(YCH4),grad(y2))*dx+y2*omegaReaction()*MolecularMass['CH4']*dx


		NewtonNS2=rho*inner(grad(u)*up,v)*dx+rho*inner(grad(up)*u,v)*dx+rhop*inner(grad(u)*u,v)*dx\
			+visco()*inner(grad(up),grad(v))*dx+inner(grad(pp),v)*dx+1/3*visco()*(v[0].dx(0)*up[0].dx(0)+v[0].dx(0)*up[1].dx(1)+v[1].dx(1)*up[0].dx(0)+v[1].dx(1)*up[1].dx(1))*dx\
#			+viscop()*inner(grad(u),grad(v))*dx-1/3*viscop()*(v[0].dx(0)*u[0].dx(0)+v[0].dx(0)*u[1].dx(1)+v[1].dx(1)*u[0].dx(0)+v[1].dx(1)*u[1].dx(1))*dx

		NewtonNS2=NewtonNS2-q*rho*rhop*inner(u,grad(rho))*dx-q*rho**2*inner(up,grad(rho))*dx-q*rho**2*inner(u,grad(rhop))*dx
		NewtonNS2=NewtonNS2-1*visco()/Pr*(4*inner(grad(rhop),grad(rho))*q+rho*inner(grad(rhop),grad(q))+rhop*inner(grad(rho),grad(q)))*dx\
#		-1/Pr*viscop()*(q*2*inner(grad(rho),grad(rho))+rho*inner(grad(rho),grad(q)))*dx
		NewtonNS2=NewtonNS2-q*(kappa-1)/kappa/patom*h0*omegaReaction()*(3*rho**2*rhop+rho**3*(\
							(Reaction['a']+Reaction['b'])*rhop/rho\
							-Reaction['E_a']/Reaction['R']*rhop/rho/(1.013e5/rho/(8314.4598/28.949))\
							+Reaction['a']*YO2p/YO2\
							+Reaction['b']*YCH4p/YCH4))*dx

	#Linearised equation : NewtonNS2
		NewtonNS2=NewtonNS2+(s*rho**3*div(up)*dx+s*3*rho**2*rhop*div(u)*dx)
		NewtonNS2=NewtonNS2-1*visco()/Pr*(4*inner(grad(rhop),grad(rho))*s+rho*inner(grad(rhop),grad(s))+rhop*inner(grad(rho),grad(s)))*dx\
#		-1/Pr*viscop()*(s*2*inner(grad(rho),grad(rho))+rho*inner(grad(rho),grad(s)))*dx
		NewtonNS2=NewtonNS2-s*(kappa-1)/kappa/patom*h0*omegaReaction()*(3*rho**2*rhop+rho**3*(\
							(Reaction['a']+Reaction['b'])*rhop/rho\
							-Reaction['E_a']/Reaction['R']*rhop/rho/(1.013e5/rho/(8314.4598/28.949))\
							+Reaction['a']*YO2p/YO2\
							+Reaction['b']*YCH4p/YCH4))*dx


		NewtonNS2=NewtonNS2+y1*rhop*inner(u,grad(YO2))*dx+y1*rho*inner(up,grad(YO2))*dx+y1*rho*inner(u,grad(YO2p))*dx
		NewtonNS2=NewtonNS2	+visco()/Sc*inner(grad(YO2p),grad(y1))*dx\
#			+viscop()*inner(grad(YO2),grad(y1))*dx
		NewtonNS2=NewtonNS2+y1*omegaReaction()*MolecularMass['O2']*2*\
									(\
							+Reaction['a']*YO2p/YO2\
							+Reaction['b']*YCH4p/YCH4\
							+(Reaction['a']+Reaction['b'])*rhop/rho\
							-Reaction['E_a']/Reaction['R']*rhop/rho/(1.013e5/rho/(8314.4598/28.949)\
							))*dx

		NewtonNS2=NewtonNS2+y2*rhop*inner(u,grad(YCH4))*dx+y2*rho*inner(up,grad(YCH4))*dx+y2*rho*inner(u,grad(YCH4p))*dx
		NewtonNS2=NewtonNS2	+visco()/Sc*inner(grad(YCH4p),grad(y2))*dx\
#			+viscop()*inner(grad(YCH4),grad(y2))*dx
		NewtonNS2=NewtonNS2+y2*omegaReaction()*MolecularMass['CH4']*\
									(\
							+Reaction['a']*YO2p/YO2\
							+Reaction['b']*YCH4p/YCH4\
							+(Reaction['a']+Reaction['b'])*rhop/rho\
							-Reaction['E_a']/Reaction['R']*rhop/rho/(1.013e5/rho/(8314.4598/28.949)\
							))*dx

		#Solve
		NewtonNS2=-NewtonNS2

		perturbation = Function(W)
		solve(NewtonNS2 == NewtonNS, perturbation, bcs)

		#Control increment
#		perturbation.vector()[W.sub(0).dofmap().dofs()]=perturbation.vector()[W.sub(0).dofmap().dofs()]*0.0
#		perturbation.vector()[W.sub(1).dofmap().dofs()]=perturbation.vector()[W.sub(1).dofmap().dofs()]*0.0
#		perturbation.vector()[W.sub(2).dofmap().dofs()]=perturbation.vector()[W.sub(2).dofmap().dofs()]*0.001
#		perturbation.vector()[W.sub(3).dofmap().dofs()]=perturbation.vector()[W.sub(3).dofmap().dofs()]*0.0
#		perturbation.vector()[W.sub(4).dofmap().dofs()]=perturbation.vector()[W.sub(4).dofmap().dofs()]*0.0
		base.vector()[:]=base.vector()[:]+perturbation.vector()[:]*0.05

		#Clip base flow
		#base.vector()[W.sub(2).dofmap().dofs()]=np.clip(base.vector()[W.sub(2).dofmap().dofs()],np.min(MF['rho'].vector()[:]),np.max(MF['rho'].vector()[:]))
		#base.vector()[W.sub(2).dofmap().dofs()]=np.clip(base.vector()[W.sub(2).dofmap().dofs()],0.1,1.16)
		base.vector()[W.sub(3).dofmap().dofs()]=np.clip(base.vector()[W.sub(3).dofmap().dofs()],1e-15,0.5)
		base.vector()[W.sub(4).dofmap().dofs()]=np.clip(base.vector()[W.sub(4).dofmap().dofs()],1e-15,0.2)
		i=i+1
		err=np.linalg.norm(perturbation.vector()[W.sub(2).dofmap().dofs()]) #error defined as the norm of density increment
		print ("Iteration "+str(i)+'. error : '+str(err))


		if i%5==0: #plot base flow and increment each 5 iterations
			fig, ax = plt.subplots()
			plt.subplot(611)
			cs=plot(base[0])
			cbar = fig.colorbar(cs)
			plt.subplot(612)
			cs=plot(base[1])
			cbar = fig.colorbar(cs)
			plt.subplot(613)
			cs=plot(base[2])
			cbar = fig.colorbar(cs)
			plt.subplot(614)
			cs=plot(base[3])
			cbar = fig.colorbar(cs)
			plt.subplot(615)
			cs=plot(base[4])
			cbar = fig.colorbar(cs)
			plt.subplot(616)
			cs=plot(base[5])
			cbar = fig.colorbar(cs)

			fig, ax = plt.subplots()
			plt.subplot(611)
			cs=plot(perturbation[0])
			cbar = fig.colorbar(cs)
			plt.subplot(612)
			cs=plot(perturbation[1])
			cbar = fig.colorbar(cs)
			plt.subplot(613)
			cs=plot(perturbation[2])
			cbar = fig.colorbar(cs)
			plt.subplot(614)
			cs=plot(perturbation[3])
			cbar = fig.colorbar(cs)
			plt.subplot(615)
			cs=plot(perturbation[4])
			cbar = fig.colorbar(cs)
			plt.subplot(616)
			cs=plot(perturbation[5])
			cbar = fig.colorbar(cs)

#Save xml
#str_E_Re='Sun1_h1_E0p755_final'
#u0,p0,rho0,yO20,yCH40=base.split()
#u_old=project(u0,VF)
#File(str_E_Re+'_u.xml') << u_old
#p_old=project(p0,RF)
#File(str_E_Re+'_p.xml') << p_old
#rho_old=project(rho0,RF)
#File(str_E_Re+'_rho.xml') << rho_old
#yO2_old=project(yO20,RF)
#File(str_E_Re+'_yO2.xml') << yO2_old
#yCH4_old=project(yCH40,RF)
#File(str_E_Re+'_yCH4.xml') << yCH4_old

#Save pvd
#u0,p0,rho0,yO20,yCH40=base.split()
#File('Sun1_h1_E0p755_final_base_flow_u.pvd') << u0
#File('Sun1_h1_E0p755_final_base_flow_rho.pvd') << rho0
#File('Sun1_h1_E0p755_final_base_flow_yO2.pvd') << yO20
#File('Sun1_h1_E0p755_final_base_flow_yCH4.pvd') << yCH40
#
#u0p,p0p,rho0p,yO20p,yCH40p=perturbation.split()
#File('Sun1_h1_E0p755_final_perturbation_flow_u.pvd') << u0p
#File('Sun1_h1_E0p755_final_perturbation_flow_rho.pvd') << rho0p
#File('Sun1_h1_E0p755_final_perturbation_flow_yO2.pvd') << yO20p
#File('Sun1_h1_E0p755_final_perturbation_flow_yCH4.pvd') << yCH40p

#plot xml
#u_new_new=Function(VF,'new_Re506_base_flow_u__6.xml')
#plot(u_new_new)
#yO2_new_new=Function(RF,'E_Re_yO2.xml')
#plot(yO2_new_new)