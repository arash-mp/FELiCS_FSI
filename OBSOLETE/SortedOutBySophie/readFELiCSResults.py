'''
Input: MeshFile, MeanFlowFile, SolutionFile (in this order)

This function reads the files created by FELiCS and writes the data to a dictionary.
Data can be plotted optionally.
'''

import hdfdict
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def readFELiCSResults(MeshFile, MeanFlowFile, SolutionFile):

	results = {}
	
	with h5py.File(MeshFile, "r") as file:
		meshDict	= {}
		cells  	= list(file.keys())[0]
		coords 	= list(file.keys())[1]
		meshDict["cells"] 	= hdf5DatasetToDict(file[cells]) 
		coordinatesDict 	= {}
		for element in list(file[coords]):
			values = hdf5DatasetToDict(file[coords])
			data  = values[element]
			coordinatesDict[element] = data
		meshDict["coordinates"] = coordinatesDict
		results["mesh"] 	= meshDict

	with h5py.File(MeanFlowFile, "r") as file:
		meanflowDict = {}
		for element in list(file['meanflow']):
			values = hdf5DatasetToDict(file['meanflow'][element])
			data = values['magnitude']		
			meanflowDict[element] = data
		results["meanflow"] = meanflowDict

	with h5py.File(SolutionFile, "r") as file:
		fluctuationDict = {}
		for element in list(file['fluctuation']['0']['pointData'].keys()):
			values = hdf5DatasetToDict(file['fluctuation']['0']['pointData'][element])
			magnitudes = values['magnitude']
			angles = values['angle']
			data = magnitudes * np.exp(angles*1j)
			fluctuationDict[element] = data
		results['fluctuation'] = fluctuationDict
		
	
	plotData = 1
	
	if plotData:
		x 		= results['mesh']['coordinates']['x']
		y 		= results['mesh']['coordinates']['y']
		triangles 	= results['mesh']['cells']['triangles']	
		triang 		= mtri.Triangulation(x, y, triangles)	

		T 		= results['meanflow']['T']
	
		plt.figure(1)
		plt.tripcolor(triang, T)
		plt.show()
	
def hdf5DatasetToDict(ds):
	keys = list(ds)
	return {keys[index]: v[()] for index, v in enumerate(list(ds.values()))}
