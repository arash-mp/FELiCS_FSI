import pdb
import numpy as np
import copy

def runCase(param, useGUI):
	'''This function runs the calculations preset in param
	Input:
		param: Parameter objects (see parameters.py), defining the case
		useGUI: Boolean, True if program is run using GUI, False if run from
		terminal directly
	'''
	from   FELiCS.ExportSolution import ExportFromFile
	import FELiCS.WeakFormulationCollection as WeakFormulationCollection
	import FELiCS.Import as Import
	import FELiCS.DefineFEMSpaces as DefineFEMSpaces
	from   FELiCS.ExportSolution import ExportGUI,ExportFromFile
	from   FELiCS.meanFlowClass import meanFlowClass
	from   FELiCS.fluctuationClass import fluctuationSolutions

	mesh=param.BCs.getMesh()
	print('Defining FEMSpaces...')
	FEMSpaces         =DefineFEMSpaces.FEMSpacesClass(
														param,
														mesh,
														)


	print('Reading InputFlow...')
	MeanFlow = meanFlowClass(param, FEMSpaces)
	MeanFlow.importDataFromFile()
	if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
		MeanFlow.exportBaseFlowAsHDF5()

	# export the mapped Meanflow:
	meanflowFilename = 'meanflow.h5'
	MeanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)

	print('Discretizing the Equations...')
	WeakFormulation = WeakFormulationCollection.WeakFormulationCollectionClass(
					param,
					FEMSpaces,
					MeanFlow
					)
	LinearAlgebraObj = WeakFormulation.DiscretizeFlow()



	if param.Case.AnalysisMode=='Modal':
		fluctSolutList = LinearAlgebraObj.\
		solveGEVP(False)
		if param.Case.CalculateAdjoint:
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
