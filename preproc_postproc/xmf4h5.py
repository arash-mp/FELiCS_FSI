import numpy as np

from os import listdir
from os.path import relpath
import sys
from fenics import Mesh
import h5py

if sys.argv[1] in ["--help","-h"]:
	print("Tool to create a xml file for FELiCS-output-h5 files for subsequent postprocessing") 
	print("Use as follows:")
	print("python xmf4h5.py --sol path/to/solution --mesh path/to/mesh")
	exit()
	

else:
	for index,arg in enumerate(sys.argv[1:]):
		if arg in ["--sol","-s"]:
			solutionPath = sys.argv[index+2]
		if arg in ["--mesh","-m"]:
			meshPath = sys.argv[index+2]

outputFolder,outputFile=solutionPath.rsplit('/', 1)
if  outputFile.find('<sol>') != -1:
	firstNamePart,secondNamePart=outputFile.split('<sol>')
	allSolutFiles=sorted([f for f in listdir(outputFolder) if f.endswith(secondNamePart) and f.startswith(firstNamePart)])
print(allSolutFiles)


if meshPath.split('.')[-1]=='xml':
	mesh = Mesh(meshPath)
	coordinates = np.array(mesh.coordinates())
	tri=mesh.cells()
	#write file in h5 format
	meshPath=meshPath.rsplit('.',1)[0]+".h5"
	hf = h5py.File(meshPath, 'w')
	h5coordinates = hf.create_group('Coordinates')
	h5coordinates.create_dataset('x',data=coordinates[:, 0])
	h5coordinates.create_dataset('y',data=coordinates[:, 1])
	h5triangulation = hf.create_group('Triangulation')
	h5triangulation.create_dataset('Triangles',data=tri.flatten())
	
# now read h5 file and forget the xml mesh file
meshh5 = h5py.File(meshPath, 'r')
tri= meshh5['Triangulation']['Triangles']

tri=np.array(tri).reshape((int(len(tri)/3),3))



for i,solutFileName in enumerate(allSolutFiles):
	solutionPath=outputFolder+ '/' +solutFileName
	solution = h5py.File(solutionPath, 'r')
	PointDataFields=list(solution['PointData'].keys())
	relMeshPath=relpath(meshPath,outputFolder)
	with open(outputFolder+'/'+solutFileName[:-2]+'xmf', 'w') as writer:
		writer.write('<?xml version="1.0" ?>\n')
		writer.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
		writer.write('  <Domain>\n')
		writer.write('     <Grid Collection="Triangle_Mesh" Name="solution-Triangle">\n')
		writer.write('     <Time Value=" '+str(i)+'" />\n')
		writer.write('        <Topology Type="Triangle" NumberOfElements="     '+str(len(tri[:,0]))+'">\n')
		writer.write('          <DataItem ItemType="Function" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'" Function="$0 - 0">\n')
		writer.write('           <DataItem Format="HDF" DataType="Int" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'">\n')
		writer.write('              '+relMeshPath+':/Triangulation/Triangles\n')
		writer.write('           </DataItem>\n')
		writer.write('          </DataItem>\n')
		writer.write('        </Topology>\n')
		writer.write('        <Geometry Type="X_Y">\n')
		for key in list(meshh5['Coordinates'].keys()):
			writer.write('           <DataItem Format="HDF" ItemType="Uniform" Precision="8" NumberType="Float" Dimensions="     '+str(len(meshh5['Coordinates'][key]))+'">\n')
			writer.write('              '+relMeshPath+':/Coordinates/'+key+'\n')
			writer.write('           </DataItem>\n')
		writer.write('        </Geometry>\n')
		for key in list(solution['PointData'].keys()):
			writer.write('       <Attribute Name="'+key+'" Center="Node" AttributeType="Scalar">\n')
			writer.write('           <DataItem Format="HDF" ItemType="Uniform" Precision="4" NumberType="Float" Dimensions="     '+str(len(solution['PointData'][key]))+'">\n')
			writer.write('              '+solutFileName+':/PointData/'+key+'\n')
			writer.write('           </DataItem>\n')
			writer.write('        </Attribute>\n')
		writer.write('     </Grid>\n')
		writer.write('  </Domain>\n')
		writer.write('</Xdmf>\n')
