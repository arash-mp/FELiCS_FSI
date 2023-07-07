#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * opermission of the copyright owner, the FLOW group at TU Berlin
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de
# */
#%%
'''
# **********************************************************************
# * This file deals with the Export of the Results data
# *
# * This file was created by Thomas L. Kaiser. Significant contributions
# * were made by
# * -
# *
# *
# ********************
'''
from sys import exit
# import tkinter as tk
# from tkinter import filedialog
#
# import matplotlib
# matplotlib.use("Agg")
# from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# import matplotlib.pyplot as plt

from dolfinx.fem import (
						Function,
)

#from fenics import Function,File,plot,cpp
from functions import *
import numpy as np
#import pyvtk
#import matplotlib.tri as tri
from colorama import Fore, Style
import scipy.io as sio

import multiprocessing as mp
from functools import partial

import pdb

def writeCSVGains(param,gains):
	''' This file writes the gains to a CSV file in the solution directory provided by the user'''
	import csv
	with open(param.Export.ExportFolder+'/gains.csv', mode='w') as writer_file:
		writer = csv.writer(writer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		n_gains,n_omega=np.shape(gains)
		head=['omega']
		for i in range(n_gains):
			head.append('gain'+str(i))
		writer.writerow(head)
		for i in range(n_omega):
			line=[param.IOResolvent.Omegas[i]]
			for j in range(n_gains):
				line.append(gains[j,i])
			writer.writerow(line)

def writeCSVSpectrum(param,spectrumDirect,spectrumAdjoint):
	''' This file writes the gains to a CSV file in the solution directory provided by the user'''
	import csv
	print (param.Export.ExportFolder+'/spectrum.csv')
	with open(param.Export.ExportFolder+'/spectrum.csv', mode='w') as writer_file:
		writer = csv.writer(writer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		head=['omega_direct_r','omega_direct_i','omega_adjoint_r','omega_adjonint_i']
		writer.writerow(head)
		for i in range(0,len(spectrumDirect)):
			line=[str(np.real(spectrumDirect[i])),str(np.imag(spectrumDirect[i])),str(np.real(spectrumAdjoint[i])),str(np.imag(spectrumAdjoint[i]))]
			writer.writerow(line)


def writeLastGitCommit(param):
	commitLabel=getLastGitCommit()
	gitFile=open(param.Export.ExportFolder+'/LastCommit.txt', mode='w')
	gitFile.write(commitLabel)


def ExportGUI(param, fluctSolutList, MeanFlow,FEMSpace, WeakFormulation):

	"""
	Function exporting the eigenvectors to a format asked for by the user


	Function arguments:
	- param:
	- fluctSolutList:

	Function returns:

	"""
	# build arrays containing the
	#A function that will be called, when the user clicks on the plot window...
	def onclick(event):
		global canvas2
		ix, iy = event.xdata, event.ydata
		if param.Case.AnalysisMode in ['Resolvent','Input-Output']:
#			xmin, xmax = ax.get_ylim()
			ix_pos=(ix-xmin)/(xmax-xmin)
			iy_pos=(np.log(iy)-np.log(ymin))/(np.log(ymax)-np.log(ymin))
			omega_xpos=(param.IOResolvent.Omegas-xmin)/(xmax-xmin)
			gains_ypos=(np.log(gains)-np.log(ymin))/(np.log(ymax)-np.log(ymin))
			n_omega = len(param.IOResolvent.Omegas)
			distance=np.ones(np.shape(gains))
			for (i,j), x in np.ndenumerate(gains_ypos):
				distance[i,j]=(omega_xpos[j]-ix_pos)**2+(x-iy_pos)**2
			## The user has clicked on the gain with the nSolution_index in n_Solution and the omega_index in omega...
			nSolut_index,omega_index=np.unravel_index(np.argmin(distance, axis=None), distance.shape)
			#nSolut_index = nSolut_index - 1
			## Read out the corresponding pair of response and forcing...
			nSolutForPlot = nSolut_index-1
			response_real.vector[:]=np.real(responses[:,nSolutForPlot,omega_index]).astype(float)
			if param.Case.AnalysisMode in ['Resolvent']:
				forcing_real.vector[:]=np.real(forcings[:,nSolutForPlot,omega_index]).astype(float)
			if param.Case.nDim==2:
				## Plot the corresponding pair of response and forcing...
				fig2 = plt.figure()
				plt.set_cmap('coolwarm')
				# For resolvent analysis we have to plot both forcing and response (2 columns), while for the IO analysis we
				# only need to plot the resposne (1 columns)
				if param.Case.AnalysisMode in ['Resolvent']:
					nPlotColumns=2
				else:
					nPlotColumns=1
				plot_index=0
				for name in param.SolutionList:
					if name == 'u':
						for component in param.VelocityComponents:
							if param.Case.AnalysisMode in ['Resolvent']:
								ax = fig2.add_subplot(nSolutionFields,nPlotColumns,plot_index*nPlotColumns+1)
								cs=plot(forcing_real.split()[0].split()[plot_index])
								cbar = fig2.colorbar(cs)
								plt.title('Forcing in '+ name+component)
							ax = fig2.add_subplot(nSolutionFields,nPlotColumns,plot_index*nPlotColumns+nPlotColumns)
							cs=plot(response_real.split()[0].split()[plot_index])
							cbar = fig2.colorbar(cs)
							plt.title('Forcing in '+ name+component)
							plot_index+=1
					else:
						if param.Case.AnalysisMode in ['Resolvent']:
							ax = fig2.add_subplot(nSolutionFields,nPlotColumns,plot_index*nPlotColumns+1)
							cs=plot(forcing_real.split()[plot_index-param.nVelocityComponents+1])
							cbar = fig2.colorbar(cs)
							plt.title('Forcing in '+ name+component)
						ax = fig2.add_subplot(nSolutionFields,nPlotColumns,plot_index*nPlotColumns+nPlotColumns)
						if len(param.SolutionList) == 1:
							#pdb.set_trace()
							#ax.text(2, 6, r'Cant plot Preview for Single Solution. PLease check the written h5. Outputs', fontsize=15)
							#responseIfOnlyOneField.vector()[:] = response_real.compute_vertex_values()
							cs = plot(responseIfOnlyOneField)
						else:
							cs=plot(response_real.split()[plot_index-param.nVelocityComponents+1])
							cbar = fig2.colorbar(cs)
						plt.title('Forcing in '+ name)
						plot_index+=1

				if len(resultWin.children)<2:
					canvas2 = FigureCanvasTkAgg(fig2, master = resultWin)
					canvas2.draw()
					tBar = NavigationToolbar2Tk( canvas2, resultWin )
					tBar.update()
					canvas2._tkcanvas.pack(side=tk.RIGHT)
				else:
					resultWin.winfo_children()[2].destroy()
					resultWin.winfo_children()[1].destroy()
					canvas2 = FigureCanvasTkAgg(fig2, master = resultWin)
					canvas2.draw()
					tBar = NavigationToolbar2Tk( canvas2, resultWin )
					tBar.update()
					canvas2._tkcanvas.pack(side=tk.RIGHT)
					canvas2.toolbar.update()
			if param.Case.AnalysisMode in ['Input-Output']:
				for omega_index in range(len(param.IOResolvent.Omegas)):

					# search the solution Objects, where the gain Number is the nSOlut_index.
					# Maybe change the property-name gainNumb to nSolutIndex?
					for fluctSolut in fluctSolutList:
						if fluctSolut.omega == param.IOResolvent.Omegas[omega_index]:
							filename = f'{param.Case.AnalysisMode}_Omega{np.round(fluctSolut.omega, 3)}_{fluctSolut.solutionKind}_gain{fluctSolut.gainNumber}.h5'
							fluctSolut.exportSolution(filename, 'o')
			else:

				# search the solution Objects, where the gain Number is the nSOlut_index.
				# Maybe change the property-name gainNumb to nSolutIndex?
				for fluctSolut in fluctSolutList:
					if fluctSolut.gainNumber == nSolut_index and fluctSolut.omega == param.IOResolvent.Omegas[omega_index]:
						filename = f'{param.Case.AnalysisMode}_Omega{np.round(fluctSolut.omega, 3)}_{fluctSolut.solutionKind}_gain{fluctSolut.gainNumber}.h5'
						fluctSolut.exportSolution(filename, 'o')
			writeLastGitCommit(param)

		if param.Case.AnalysisMode == 'Modal':
			ix_pos=(ix-xmin)/(xmax-xmin)
			iy_pos=(iy-ymin)/(ymax-ymin)
			# normalizes the plotted eigenavalues to window size to suit the format of mouse click
			EValDirect_xpos=(np.real(eValDirect)-xmin)/(xmax-xmin)
			EValDirect_ypos=(np.imag(eValDirect)-ymin)/(ymax-ymin)
			EValAdj_xpos=(np.real(eValAdjoint)-xmin)/(xmax-xmin)
			EValAdj_ypos=(np.imag(eValAdjoint)-ymin)/(ymax-ymin)
			# distance from mouse click to direct eval and adjoint eval
			distanceDir=np.ones(np.shape(eValDirect))
			distanceAdj=np.ones(np.shape(eValDirect))
			for i, x in np.ndenumerate(distanceDir):
				distanceDir[i]=(EValDirect_xpos[i]-ix_pos)**2+(EValDirect_ypos[i]-iy_pos)**2
				distanceAdj[i]=(EValAdj_xpos[i]-ix_pos)**2+(EValAdj_ypos[i]-iy_pos)**2

			# distingueshes wether a direct or adjoint eval was close to mouse click
			# finds out first index -> determines the corresponding Eval
			# substracts the conjugate from the other eval-array and finds minimum
			# which corresponds to the twin of the first eval found
			if np.min(distanceDir)<np.min(distanceAdj):
				EValDirIndex=np.unravel_index(np.argmin(distanceDir, axis=None), distanceDir.shape)[0]
				EValAdjIndex = np.argmin( np.abs( eValAdjoint\
							  - np.conj(eValDirect[EValDirIndex])) )

			elif np.min(distanceDir)>=np.min(distanceAdj):
				EValAdjIndex=np.unravel_index(np.argmin(distanceAdj, axis=None), distanceAdj.shape)[0]
				EValDirIndex = np.argmin( np.abs( eValDirect\
							  - np.conj(eValAdjoint[EValAdjIndex])) )



    		## Write out the corresponding Direct and Adjoint mode
			ModeDirect.vector[:]=np.real(eVecDirect[:,EValDirIndex].flatten()).astype(float)
			ModeAdjoint.vector[:]=np.real(eVecAdjoint[:,EValAdjIndex].flatten()).astype(float)
    		#Plot the modes chosen by the user...
			fig2 = plt.figure()
			plt.set_cmap('coolwarm')
			plot_index=0
			for name in param.SolutionList:
				if name =='u':
					for component in param.VelocityComponents:
						ax = fig2.add_subplot(nSolutionFields,2,plot_index*2+1)
						cs=plot(ModeDirect.split()[0].split()[plot_index])
						cbar = fig2.colorbar(cs)
						plt.title('Direct Mode in '+ name+component)
						ax = fig2.add_subplot(nSolutionFields,2,plot_index*2+2)
						cs=plot(ModeAdjoint.split()[0].split()[plot_index])
						cbar = fig2.colorbar(cs)
						plt.title('Adjoint Mode in '+ name+component)
						plot_index += 1
				else:
					ax = fig2.add_subplot(nSolutionFields,2,plot_index*2+1)

					cs=plot(ModeDirect.split()[plot_index-param.nVelocityComponents+1])
					cbar = fig2.colorbar(cs)
					plt.title('Direct Mode in '+ name+component)
					ax = fig2.add_subplot(nSolutionFields,2,plot_index*2+2)
					cs=plot(ModeAdjoint.split()[plot_index-param.nVelocityComponents+1])
					cbar = fig2.colorbar(cs)
					plt.title('Adjoint Mode in '+ name+component)
					plot_index += 1

			if len(resultWin.children)<2:
				canvas2 = FigureCanvasTkAgg(fig2, master = resultWin)
				canvas2.draw()
				tBar = NavigationToolbar2Tk( canvas2, resultWin )
				tBar.update()
				canvas2._tkcanvas.pack(side=tk.RIGHT)
			else:
				resultWin.winfo_children()[2].destroy()
				resultWin.winfo_children()[1].destroy()
				canvas2 = FigureCanvasTkAgg(fig2, master = resultWin)
				canvas2.draw()
				tBar = NavigationToolbar2Tk( canvas2, resultWin )
				tBar.update()
				canvas2._tkcanvas.pack(side=tk.RIGHT)
				canvas2.toolbar.update()

			for fluctSolut in fluctSolutList:

				if eValDirect[EValDirIndex] == fluctSolut.omega or eValAdjoint[EValAdjIndex] == fluctSolut.omega:

					filename = f'{param.Case.AnalysisMode}Solution_Omega_{fluctSolut.solutionKind}_{np.round(fluctSolut.omega, 3)}.h5'
					fluctSolut.exportSolution(filename, 'o')
			writeLastGitCommit(param)


#		resultWin.destroy()
	def reset_ix_glob():
		global ix_glob, iy_glob
		del ix_glob,iy_glob
	#index=0
	# Get the Names of the Solution Fields (ux, uy etc...)
	#SolutionStrings=getSolutionInfo(param)
	SolutionStrings = param.SolutionList
	# Number of Solutions to write
	nSolutionFields = len(SolutionStrings)+param.nVelocityComponents-1
	# Get the mesh data...
	mesh=FEMSpace.P2.mesh
	xy = mesh.coordinates()
	mesh_cells=mesh.cells()
#	plt.ioff()
	#If user chose resolvent analysis...
	if param.Case.AnalysisMode in ['Resolvent','Input-Output']:

		# get the gains from the list of fluctuation Solutions:
		gains = np.zeros((param.Numerics.nSolut, len(param.IOResolvent.Omegas)), 'complex')
		# get nDOFS:
		nDOF = fluctSolutList[0].solutVector.shape[0]
		responses = []
		forcings = []
		responses = np.zeros((nDOF,param.Numerics.nSolut, len(param.IOResolvent.Omegas)),'complex')
		forcings = np.zeros((nDOF,param.Numerics.nSolut, len(param.IOResolvent.Omegas)),'complex')

		for i, fluctSolut in enumerate(fluctSolutList):
			indexOmega = param.IOResolvent.Omegas.index(fluctSolut.omega)
			if fluctSolut.solutionKind == 'Response':
				gains[fluctSolut.gainNumber, indexOmega] = fluctSolut.gainValue
				responses[:, fluctSolut.gainNumber, indexOmega] = fluctSolut.solutVector
			else:
				forcings[:, fluctSolut.gainNumber, indexOmega] = fluctSolut.solutVector



		writeCSVGains(param,gains)
		#Define Function Spaces for the solutions(response and forcing)
		response_real = Function(FEMSpace.VMixed)

		responseIfOnlyOneField = Function(FEMSpace.P1)

		response_imag = Function(FEMSpace.VMixed)
		response_abs = Function(FEMSpace.VMixed)
		response_ang = Function(FEMSpace.VMixed)
		forcing_real = Function(FEMSpace.VMixed)
		forcing_imag = Function(FEMSpace.VMixed)
		forcing_abs = Function(FEMSpace.VMixed)
		forcing_ang = Function(FEMSpace.VMixed)
		# While InPlotLoop ==true stay in output loop (So far the loop cannot be ended. This should be changed some day)
		resultWin = tk.Tk()
		resultWin.title('Results Plot')

		### Plot direct and adjoint spectrum
		fig1 = plt.figure()
		ax = fig1.add_subplot(1,1,1)
		ax.plot(param.IOResolvent.Omegas,gains.transpose(),'x')
		plt.xlabel(r'$\omega_r$')
		plt.ylabel(r'gain')
		plt.title('left: 2D-single, left double: 2D-video, right: 3D-single, right double: 3D-video')

#		plt.legend()
		ymin, ymax = ax.get_ylim()
		ax.set_ylim(0.1,ymax)
		ymin, ymax = ax.get_ylim()
		ax.autoscale(enable=True, axis='both', tight=None)
		xmin, xmax = ax.get_xlim()
		ax.set_yscale('log')
		plt.grid()

		canvas = FigureCanvasTkAgg(fig1, master = resultWin)
		canvas.get_tk_widget().pack(side=tk.LEFT)
		canvas.mpl_connect('button_press_event', onclick)
		resultWin.mainloop()



	#If user chose resolvent analysis...
	elif param.Case.AnalysisMode == 'Modal':

		# extract the eigenvalues of the direct- and adjoint Problem from the
		# fluctuation solution List:
		eValDirect = []
		eValAdjoint = []
		eVecDirect = []
		eVecAdjoint = []

		nDof = fluctSolutList[0].solutVector.shape[0]
		eVecDirect = np.zeros((nDof, param.Numerics.nSolut*len(param.Numerics.EigenValueGuess)), 'complex')
		eVecAdjoint = np.zeros((nDof, param.Numerics.nSolut*len(param.Numerics.EigenValueGuess)), 'complex')
		directIndex = 0
		adjointIndex = 0
		for fluctSolut in fluctSolutList:
			if fluctSolut.solutionKind == 'Direct':

				eValDirect.append(fluctSolut.omega)
				eVecDirect[:, directIndex] = fluctSolut.solutVector
				directIndex += 1
			else:
				eValAdjoint.append(fluctSolut.omega)
				eVecAdjoint[:, adjointIndex] = fluctSolut.solutVector
				adjointIndex += 1

		# convert the list of eVals to a np-array:
		eValDirect = np.array(eValDirect)
		eValAdjoint = np.array(eValAdjoint)

		writeCSVSpectrum(param,eValDirect,eValAdjoint)
		writeLastGitCommit(param)
		## Define Function spaces for the direct mode and the adjoint mode
		ModeDirect = Function(FEMSpace.VMixed)
		ModeAdjoint = Function(FEMSpace.VMixed)
		## While InPlotLoop == True remain in export loop (As in the resolvent analysis a way to exit thios loop remains to be implemented...)
		resultWin = tk.Tk()
		resultWin.title('Results Plot')

		### Plot direct and adjoint spectrum
		fig1 = plt.figure()
		ax = fig1.add_subplot(1,1,1)
		ax.plot(np.real(eValDirect[:]),np.imag(eValDirect[:]),'x', label='Direct')
		ax.plot(np.real(eValAdjoint[:]),np.imag(eValAdjoint[:]),'o', label='Adjoint')
		# Here the function for getting the user click position is called

		plt.ylabel(r'$\omega_i$')
		plt.xlabel(r'$\omega_r$')
		plt.title('left: 2D-single, left double: 2D-video, right: 3D-single, right double: 3D-video')
		plt.legend()
		ymin, ymax = ax.get_ylim()
		xmin, xmax = ax.get_xlim()

		canvas = FigureCanvasTkAgg(fig1, master = resultWin)
		canvas.get_tk_widget().pack(side=tk.LEFT)
		canvas.mpl_connect('button_press_event', onclick)
		resultWin.mainloop()

def ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow):
	# automated export from here on-Felics call wihtout gui --> no export-window: all frequencies and the first 2 leading modes
	xy = FEMSpaces.P2.mesh.coordinates()
	mesh_cells=FEMSpaces.P2.mesh.cells()
	# in the case of resolvent analysis: expport all frequencies and the first two leading modes

	if param.Case.AnalysisMode in ['Resolvent','Input-Output']:


		# get the gains from the list of fluctuation Solutions:
		gains = np.zeros((param.Numerics.nSolut, len(param.IOResolvent.Omegas)))

		ct = 1
		for fluctSolut in fluctSolutList:
			print(f"-- Saving output file {ct:0.0f} / {len(fluctSolutList):0.0f}")
			if fluctSolut.solutionKind == 'Response':
				indexOmega = param.IOResolvent.Omegas.index(fluctSolut.omega)
				gains[fluctSolut.gainNumber, indexOmega] = fluctSolut.gainValue

			fileName = f'{param.Case.AnalysisMode}_Omega{np.round(fluctSolut.omega, 3)}_{fluctSolut.solutionKind}_gain{fluctSolut.gainNumber}.h5'
			fluctSolut.exportSolution(fileName, 'o')
			ct += 1
		writeCSVGains(param,gains)
		writeLastGitCommit(param)


	# in case of modal analysis: export the whole spectrum and every mode
	elif param.Case.AnalysisMode == 'Modal':

		eValDirect = []
		eValAdjoint = []

		for fluctSolut in fluctSolutList:
			if fluctSolut.solutionKind == 'Direct':
				eValDirect.append(fluctSolut.omega)

			else:
				eValAdjoint.append(fluctSolut.omega)

		# convert the list of eVals to a np-array:
		eValDirect = np.array(eValDirect)
		eValAdjoint = np.array(eValAdjoint)

		idxDirect  = np.argmax(np.imag(eValDirect))
		idxAdjoint = np.argmax(np.imag(eValAdjoint))

		writeCSVSpectrum(param,eValDirect,eValAdjoint)

		for fluctSolut in fluctSolutList:

			fileName = f'{param.Case.AnalysisMode}Solution_Omega_{fluctSolut.solutionKind}_{np.round(fluctSolut.omega, 3)}.h5'
			fluctSolut.exportSolution(fileName, 'o')
		writeLastGitCommit(param)
