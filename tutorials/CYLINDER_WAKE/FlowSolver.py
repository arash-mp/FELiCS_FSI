#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:26:22 2019

@author: tkaiser
"""

from dolfin import *
import pyvtk
import numpy as np
#import pickle
from numpy.random import seed, random
import matplotlib.pyplot as plt
import meshio
import os

Re=47
D=1.0
U=1.0
meshfile="CylinderWakeMesh.msh"
nu=U*D/Re
tol=1.0e-10
print("LoadingMesh...")

#mesh=Mesh("CylinderWakeMesh.xml")
gmsh = meshio.read(meshfile)
os.system('dolfin-convert '+meshfile+' '+meshfile[:-4]+'.xml')
mesh = Mesh(meshfile[:-4]+'.xml')
lineIDs = np.unique(gmsh.cell_data['line']['gmsh:physical'])


## Define Boundary conditions
def symmetry(x, on_boundary):
    return x[1] < 0.000001 and on_boundary 
def inlet(x, on_boundary):
    return x[0] < (-99.999999999) and on_boundary 
def outlet(x, on_boundary):
    return x[0] > (199.999999999) and on_boundary 
def wall(x, on_boundary):
    return ((x[0]**2+x[1]**2)**0.5)<0.5000000001 and on_boundary 
def outlet(x, on_boundary):
    return x[0] > 199.999999999 and on_boundary 
def misc(x, on_boundary):
    return x[1] > 24.999999999 and on_boundary 

## Define Finite Element Spaces
P1 = FiniteElement('P', triangle, 1)
P2 = FiniteElement('P', triangle, 2)
Mixed = MixedElement([P2, P2, P1])
VMixed = FunctionSpace(mesh, Mixed)
VP2= FunctionSpace(mesh, P2)
VP1= FunctionSpace(mesh, P1)
FunctionSpaceVectorVelocity = VectorFunctionSpace(mesh, "CG", 2)
ux0ptemp= Function(VP2)

## Define Variables
vx0,vr0,q0 = TestFunctions(VMixed)
u = Function(VMixed)
u_n = Function(VMixed)
ux0,ur0,p0=split(u)
ux0p,ur0p,p0p=split(u_n)

##Assign start solution
init = Expression('1-exp(-(sqrt(x[0]*x[0]+x[1]*x[1])-.5)/.5)',degree=2)
ux0ptemp.interpolate(init)
assigner = FunctionAssigner(VMixed.sub(0), VP2)
assigner.assign(u_n.sub(0), ux0ptemp)



NewtonNS=  ((ur0*ux0p.dx(1) + ur0p*ux0.dx(1) + ux0*ux0p.dx(0) + ux0p*ux0.dx(0))*vx0 \
              + (ur0*ur0p.dx(1) + ur0p*ur0.dx(1) + ux0*ur0p.dx(0) + ux0p*ur0.dx(0))*vr0 \
              + nu * ( ux0.dx(0)*vx0.dx(0) + ux0.dx(1)*vx0.dx(1)+ \
                       ur0.dx(0)*vr0.dx(0) + ur0.dx(1)*vr0.dx(1) \
                       ) \
              - p0*( vx0.dx(0) + vr0.dx(1) ) \
              - q0*( ur0.dx(1) + ux0.dx(0) ) \
              - p0*q0*(0.000001))  *dx \
              +((ur0p*ux0p.dx(1) + ux0p*ux0p.dx(0))*vx0 \
                + ((ur0p*ur0p.dx(1) + ux0p*ur0p.dx(0)))*vr0 \
                + nu * ( ux0p.dx(0)*vx0.dx(0) + ux0p.dx(1)*vx0.dx(1)+ \
                         ur0p.dx(0)*vr0.dx(0) + ur0p.dx(1)*vr0.dx(1) \
                         ) \
                - p0p*( vx0.dx(0) + vr0.dx(1) ) \
                - q0*( ur0p.dx(1) + ux0p.dx(0) ) \
                - p0p*q0*(0.000001)) *dx
              

## Impose BCs
boundaries = MeshFunction('size_t', mesh, meshfile.split('.')[0] + '_facet_region.xml')
bcInlet1=DirichletBC(VMixed.sub(0), 0, boundaries, 1)
bcInlet2=DirichletBC(VMixed.sub(1), 0, boundaries, 1)

bcSym=DirichletBC(VMixed.sub(1), 0, boundaries, 2)

bcCyl1=DirichletBC(VMixed.sub(0), 0, boundaries, 5)
bcCyl2=DirichletBC(VMixed.sub(1), 0, boundaries, 5)
#bcInlet1=DirichletBC(VMixed.sub(0), 0.0, inlet)
#bcInlet2=DirichletBC(VMixed.sub(1), 0.0, inlet)
#bcCyl1=DirichletBC(VMixed.sub(0), 0.0, wall)
#bcCyl2=DirichletBC(VMixed.sub(1), 0.0, wall)
#bcsym=DirichletBC(VMixed.sub(1), 0.0, symmetry)
bcs=[bcInlet1,bcInlet2,bcCyl1,bcCyl2,bcSym]
    
i=0
while (i<9):
    solve(NewtonNS == 0, u, bcs)
    u_n.vector()[:]=u_n.vector()[:]+u.vector()[:]
    i=i+1
    print ("Iteration "+str(i))


action="write"
temp=u_n.split()[0].vector()[:]
u_n.split()[0].vector()[:]=u_n.split()[1].vector()[:]
u_n.split()[1].vector()[:]=temp
del temp

indexx=VMixed.sub(0).dofmap().dofs()
indexy=VMixed.sub(1).dofmap().dofs()

vel1=Function(VP2)
vel2=Function(VP2)
press=Function(VP1)
print(u_n.function_space())
assigner = FunctionAssigner([VP2,VP2,VP1], u_n.function_space())

assigner.assign([vel1,vel2,press],u_n)
u.vector()[:] = random(VMixed.dim())
vel = project(Expression(("u1", "u2"),\
	u1=vel1,\
	u2=vel2,\
	degree=2), FunctionSpaceVectorVelocity)
hdf5file = HDF5File(mesh.mpi_comm(),"cylinder_base_flow.hdf5", 'w')
#hdf5file.write(vel1, "ux")
#hdf5file.write(vel2, "uy")
hdf5file.write(vel, "u")
hdf5file.write(press, "p")
del hdf5file


fig, ax = plt.subplots()
plt.subplot(311)
cs=plot(vel1)
cbar = fig.colorbar(cs)
plt.subplot(312)
cs=plot(vel2)
cbar = fig.colorbar(cs)
plt.subplot(313)
cs=plot(press)
cbar = fig.colorbar(cs)
plt.show()

