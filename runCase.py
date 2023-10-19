from ExportSolution import ExportFromFile
import WeakFormulationCollection
import pdb

def runCase(param, useGUI):
	'''This function runs the calculations preset in param
	Input:
		param: Parameter objects (see parameters.py), defining the case
		useGUI: Boolean, True if program is run using GUI, False if run from
		terminal directly
	'''
	import Import
	import DefineFEMSpaces
	from ExportSolution import ExportGUI,ExportFromFile
	from meanFlowClass import meanFlowClass
	import global_variables as glob
	import numpy as np
	from fluctuationClass import fluctuationSolutions
	import copy

	mesh=param.BCs.getMesh()
	print('Defining FEMSpaces...')
	FEMSpaces         =DefineFEMSpaces.FEMSpacesClass(
														param,
														mesh,
														)


	print('Reading InputFlow...')
	MeanFlow = meanFlowClass(param, FEMSpaces)

	print('Discretizing the Equations...')
	WeakFormulation = WeakFormulationCollection.WeakFormulationCollectionClass(
					param,
					FEMSpaces,
					MeanFlow
					)
	LinearAlgebraObj = WeakFormulation.DiscretizeFlow()


	# export the mapped Meanflow:
	meanflowFilename = 'meanflow.h5'
	MeanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)

	if param.Case.AnalysisMode=='Modal':
		fluctSolutList = LinearAlgebraObj.\
		solveGEVP(False)
		if param.Case.AdjointFlag is True:
			fluctSolutObjListAdjointGEVP = LinearAlgebraObj.\
			solveGEVP(True)
			fluctSolutList.extend(fluctSolutObjListAdjointGEVP)

	elif param.Case.AnalysisMode=='Resolvent':
		fluctSolutList = LinearAlgebraObj.solveResolvent(
																WeakFormulation)

	elif param.Case.AnalysisMode=='Input-Output':

		fluctSolutList = LinearAlgebraObj.solveInputOutput()

	if useGUI:
		ExportGUI(param, fluctSolutList, MeanFlow,FEMSpaces, WeakFormulation)
	else:
		ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow)
