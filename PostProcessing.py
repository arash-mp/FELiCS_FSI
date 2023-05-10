#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 12:15:31 2020

@author: cwang
"""

def TransferFunctionPoint(param,FEMSpace,mean,index1,index2, response_abs, response_arg):
	# param.FlowMode in ['LowMachEnthalpy']
	#input: response_abs, response_arg, MeanFlow(for normalization),
	#       probe_names: variable names, probe_points: probe positions
	MeanFlow=mean.MeanFlowDict
	probe_names = ['ux','rho']
	probe_points = [(2e-3,0), (4e-3,0), (6e-3,0), (2e-3,1e-3), (4e-3,1e-3), (6e-3,1e-3)]

	#find index of probe names in function space
	probe_names_index=[]
	for name in probe_names:
		if name[0] == 'u':
			probe_names_index.append(param.VelocityComponents.index(name[1:]))
		else:
			probe_names_index.append(param.SolutionList.index(name)+len(param.VelocityComponents)-1)

	import os
	from numpy import pi
	def writeTF2CSV(p_name,abs_li,arg_li, index2):
		# Max: if file already exists, open in append mode, else open in write mode
		csv_name = param.Export.ExportFolder+"/TF/"+p_name+".csv"
		if os.path.isfile(csv_name):
			csv_file = open(csv_name, "a")
		else:
			os.makedirs(param.Export.ExportFolder+"/TF/",exist_ok=True)
			csv_file = open(csv_name, "w")
		csv_file.write("%s, " % param.IOResolvent.Omegas[index2])
		for value_abs in abs_li:
			csv_file.write("%s, " % value_abs)
		for value_arg in arg_li:
			csv_file.write("%s" % value_arg)
			if not value_arg == arg_li[-1]:
				csv_file.write(', ')

		csv_file.write("\n")
		csv_file.close()

	for p_ind in probe_names_index:
		abs_li = [] #temp abs list
		arg_li = [] #temp arg list
		p_name = probe_names[probe_names_index.index(p_ind)]
		if p_name[0] == 'u':
			p_MF = MeanFlow['u'][p_ind]
		else:
			p_MF = MeanFlow[p_name]
		for p_pos in probe_points:
			abs_li.append(response_abs[p_ind](p_pos)/p_MF(p_pos))
			arg_li.append(-response_arg[p_ind](p_pos)/pi*180)
		writeTF2CSV(p_name,abs_li,arg_li,index2)

def FTF(param,mean,fluc_real,fluc_imag):
	from fenics import Measure,inner,FacetNormal,dx,assemble,Constant
	import numpy as np
	MeanFlowDict=mean.MeanFlowDict
	mesh=MeanFlowDict[list(MeanFlowDict.keys())[0]].function_space().mesh()
	boundaries = param.BCs.getBoundaries()
	ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
	forcingBoundaryIndices=param.IOResolvent.ForcingBoundaryIndices
	n_BC=FacetNormal(mesh)
	S_inlet=0
	u_real=0
	u_imag=0
	u_mean=0
	for index in forcingBoundaryIndices:
		S_inlet += assemble(Constant(1) *ds(index))
	for index in forcingBoundaryIndices:
		u_real += assemble(-inner(fluc_real.u,n_BC)*ds(index))/S_inlet
		u_imag += assemble(-inner(fluc_imag.u,n_BC)*ds(index))/S_inlet
		u_mean += assemble(-inner(mean.u,n_BC)*ds(index))/S_inlet
	u_abs = np.sqrt(u_real**2+u_imag**2)
	u_ang = np.angle(np.complex(u_imag+u_imag))

	dQ_real_int = assemble(fluc_real.Q * dx)
	dQ_imag_int = assemble(fluc_imag.Q * dx)

	dQ_mean_int = assemble(mean.Q * dx)

	print(dQ_real_int)
	print(dQ_imag_int)
	dQ_abs = np.sqrt(dQ_real_int**2+dQ_imag_int**2)
	dQ_ang = np.angle(np.complex(dQ_real_int,dQ_imag_int))

	FTF_gain = (dQ_abs / dQ_mean_int) / (u_abs/u_mean)
	print(dQ_abs)
	print(dQ_mean_int)
	print(u_abs)
	print(u_mean)
	FTF_phase = np.angle(np.complex(dQ_real_int, dQ_imag_int)/np.complex(u_real, u_imag))*180/np.pi
	print('test')
	print('FTF gain: '+str(FTF_gain)+', FTF phase: '+str(FTF_phase))

def FTFAndHeatRelease(param,FEMSpace,MeanFlow,index1,index2,response_real, response_imag):
	#input: response_real, response_imag, MeanFlow(for normalization),
	#       probe_names: variable names
	#Currently implemented for Global
	hrDict={} #dict for field for heat release
	if param.Case.Mixture.ReactionMechanism['type']=='WestbrookDryer_Max':
		probe_names = ['rho','CH4']

		#find index of probe names in function space
		probe_names_index=[]
		for name in probe_names:
			if name[0] == 'u':
				probe_names_index.append(param.VelocityComponents.index(name[1:]))
			else:
				probe_names_index.append(param.SolutionList.index(name)+len(param.VelocityComponents)-1)

		from fenics import project,Function,assemble,dx
		import numpy as np
		rho_r = project(response_real.sub(probe_names_index[probe_names.index('rho')]-1), FEMSpace.P2)
		CH4_r = project(response_real.sub(probe_names_index[probe_names.index('CH4')]-1), FEMSpace.P2)

		rho_i = project(response_imag.sub(probe_names_index[probe_names.index('rho')]-1), FEMSpace.P2)
		CH4_i = project(response_imag.sub(probe_names_index[probe_names.index('CH4')]-1), FEMSpace.P2)

		from Reactions.GlobalReaction import GlobalReaction
		ReactionModelName="WestbrookDryer_Max" #to be implemented in param
		Reaction=GlobalReaction(param.Case.Mixture.ReactionMechanism)
		MeanFieldReaction=Reaction.computeMeanField(MeanFlow,FEMSpace.P2)
		HeatRelease_real=Reaction.postHeatRelease(MeanFlow, rho_r, CH4_r, FEMSpace.P2)
		HeatRelease_imag=Reaction.postHeatRelease(MeanFlow, rho_i, CH4_i, FEMSpace.P2)
		uamplitude=1.0 #hardcored. To change with Line 127
	elif param.Case.Mixture.ReactionMechanism['type']=='TwoStep':
		probe_names = ['rho','CH4','O2','CO','CO2']



		#find index of probe names in function space
		probe_names_index=[]
		for name in probe_names:
			if name[0] == 'u':
				probe_names_index.append(param.VelocityComponents.index(name[1:]))
			else:
				probe_names_index.append(param.SolutionList.index(name)+len(param.VelocityComponents)-1)

		from fenics import project,Function,assemble,dx
		import numpy as np
		rho_r = project(response_real.sub(probe_names_index[probe_names.index('rho')]-1), FEMSpace.P2)
		CH4_r = project(response_real.sub(probe_names_index[probe_names.index('CH4')]-1), FEMSpace.P2)
		O2_r = project(response_real.sub(probe_names_index[probe_names.index('O2')]-1), FEMSpace.P2)
		CO_r = project(response_real.sub(probe_names_index[probe_names.index('CO')]-1), FEMSpace.P2)
		CO2_r = project(response_real.sub(probe_names_index[probe_names.index('CO2')]-1), FEMSpace.P2)

		rho_i = project(response_imag.sub(probe_names_index[probe_names.index('rho')]-1), FEMSpace.P2)
		CH4_i = project(response_imag.sub(probe_names_index[probe_names.index('CH4')]-1), FEMSpace.P2)
		O2_i = project(response_imag.sub(probe_names_index[probe_names.index('O2')]-1), FEMSpace.P2)
		CO_i = project(response_imag.sub(probe_names_index[probe_names.index('CO')]-1), FEMSpace.P2)
		CO2_i = project(response_imag.sub(probe_names_index[probe_names.index('CO2')]-1), FEMSpace.P2)

		from Reactions.TwoStepReaction import TwoStepReaction
		ReactionModelName="BFER" #to be implemented in param
		Reaction=TwoStepReaction(ReactionModelName)
		Reaction.computeMeanField(MeanFlow.MeanFlowDict,FEMSpace.P2)
		HeatRelease_real=Reaction.postHeatRelaese(MeanFlow.MeanFlowDict, rho_r, CH4_r, O2_r, CO_r, CO2_r, FEMSpace.P2)
		HeatRelease_imag=Reaction.postHeatRelaese(MeanFlow.MeanFlowDict, rho_i, CH4_i, O2_i, CO_i, CO2_i, FEMSpace.P2)
		uamplitude=0.01 #hardcored. To change with Line 127
	else:
		print('ChemitryModel not understood. FTF and hr is not computed.')
		return 0

	#get velocity ratio
	#to discuss with Thomas
# 	err=1e-6
# 	inl=-0.03+err #inlet x position
# 	inh=0.006-err #inlet y width
# 	divi=100
# 	u_base=0
# 	u_fluc=0
# 	for ii in range(divi):
# 		u_base_tmp=MeanFlow.u[0]((inl,ii*inh/divi))
# 		u_fluc_r_tmp=response_real.sub(0)((inl,ii*inh/divi))[0]
# 		u_fluc_i_tmp=response_imag.sub(0)((inl,ii*inh/divi))[0]
# 		u_fluc_abs_tmp=np.sqrt(u_fluc_r_tmp**2+u_fluc_i_tmp**2)
# 		u_base+=u_base_tmp
# 		u_fluc+=u_fluc_abs_tmp
# 	print('uratio')
# 	print(u_fluc/u_base)
# 	print('u_fluc')
# 	print(u_fluc/divi)
# 	print('u_base')
# 	print(u_base/divi)


	HeatRelease_abs = Function(FEMSpace.P2)
	HeatRelease_abs.vector()[:] = np.sqrt(HeatRelease_real.vector()[:]**2+HeatRelease_imag.vector()[:]**2)
	HeatRelease_ang = Function(FEMSpace.P2)
	HeatRelease_ang.vector()[:] = np.angle(HeatRelease_real.vector()[:]+1j* HeatRelease_imag.vector()[:])

	hrDict['HeatRelease_real']=HeatRelease_real
	hrDict['HeatRelease_imag']=HeatRelease_imag
	hrDict['HeatRelease_abs']=HeatRelease_abs
	hrDict['HeatRelease_ang']=HeatRelease_ang

	#compute FFT
	if 'dQ' in MeanFlow.MeanFlowDict.keys():
		global_heat_s=assemble(MeanFlow.dQ*dx)
	else:
		MeanFieldReaction.vector()[:]=np.abs(MeanFieldReaction.vector()[:]*Reaction.h0)
		global_heat_s=MeanFieldReaction
	global_heat_f_real=assemble(HeatRelease_real*dx)/global_heat_s/uamplitude
	global_heat_f_imag=assemble(HeatRelease_imag*dx)/global_heat_s/uamplitude
	global_heat_f_abs=np.sqrt(global_heat_f_real**2+global_heat_f_imag**2)
	global_heat_f_ang=-np.angle(global_heat_f_real+1j*global_heat_f_imag, deg=True)

	import os
	p_name='FTF'
	csv_name = param.Export.ExportFolder+"/TF/"+p_name+".csv"
	if os.path.isfile(csv_name):
		csv_file = open(csv_name, "a")
	else:
		os.makedirs(param.Export.ExportFolder+"/TF/",exist_ok=True)
		csv_file = open(csv_name, "w")
	csv_file.write("%s, " % param.IOResolvent.Omegas[index2])
	csv_file.write("%s, " % global_heat_f_real)
	csv_file.write("%s," % global_heat_f_imag)
	csv_file.write("%s," % global_heat_f_abs)
	csv_file.write("%s" % global_heat_f_ang)
	csv_file.write("\n")
	csv_file.close()


	return hrDict
