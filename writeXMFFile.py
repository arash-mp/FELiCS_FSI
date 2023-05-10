import h5py
from os.path import relpath
def writeXMFFile(meshFile,fieldsFile):

	if '/' in fieldsFile[0]:
		outputFolder,outputFile=fieldsFile[0].rsplit('/', 1)
	else:
		outputFolder='.'
		outputFile=fieldsFile

	meshh5 = h5py.File(meshFile, 'r')
	tri= meshh5['cells']['triangles']
		

	relMeshPath=relpath(meshFile,outputFolder)
	for i,solutFileName in enumerate(fieldsFile):
		solutionPath = solutFileName
		solution = h5py.File(solutionPath, 'r')
		TimeStamps=list(solution['fluctuation'].keys())
		for i_time,time in enumerate(TimeStamps):
			writer = open(solutFileName[:-3]+'_'+str(i_time)+'.xmf', 'w')
			writer.write('<?xml version="1.0" ?>\n')
			writer.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
			writer.write('  <Domain>\n')
			writer.write('     <Grid Collection="Triangle_Mesh'+str(time)+'" Name="solution-Triangle'+str(time)+'">\n')
			writer.write('     <Time Value=" '+str(time)+'" />\n')
			writer.write('        <Topology Type="Triangle" NumberOfElements="     '+str(len(tri[:,0]))+'">\n')
			writer.write('          <DataItem ItemType="Function" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'" Function="$0 - 0">\n')
			writer.write('           <DataItem Format="HDF" DataType="Int" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'">\n')
			writer.write('              '+relMeshPath+':/cells/triangles\n')
			writer.write('           </DataItem>\n')
			writer.write('          </DataItem>\n')
			writer.write('        </Topology>\n')
			writer.write('        <Geometry Type="X_Y">\n')
			for key in list(meshh5['coordinates'].keys()):
				writer.write('           <DataItem Format="HDF" ItemType="Uniform" Precision="8" NumberType="Float" Dimensions="     '+str(len(meshh5['coordinates'][key]))+'">\n')
				writer.write('              '+relMeshPath+':/coordinates/'+key+'\n')
				writer.write('           </DataItem>\n')
			writer.write('        </Geometry>\n')
			for key in list(solution['fluctuation'][time]['pointData'].keys()):
				key_loc=key
				if 'angle' in list(solution['fluctuation'][time]['pointData'][key].keys()):
					suffix = '_magnitude'
				else:
					suffix = ''
				
				#Write magnitude
				writer.write('       <Attribute Name="'+key_loc+suffix+'" Center="Node" AttributeType="Scalar">\n')
				writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution['fluctuation'][time]['pointData'][key]['angle']))+'" Function="$0">\n')
				writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['magnitude']))+'">\n')
				writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/magnitude\n')
				writer.write('               </DataItem>\n')
				writer.write('           </DataItem>\n')
				writer.write('        </Attribute>\n')
				if 'angle' in list(solution['fluctuation'][time]['pointData'][key].keys()):
					#Write angle 
					writer.write('       <Attribute Name="'+key_loc+'_angle" Center="Node" AttributeType="Scalar">\n')
					writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution['fluctuation'][time]['pointData'][key]['angle']))+'" Function="$0">\n')
					writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['angle']))+'">\n')
					writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/angle\n')
					writer.write('               </DataItem>\n')
					writer.write('           </DataItem>\n')
					writer.write('        </Attribute>\n')
					#Write real part
					writer.write('       <Attribute Name="'+key_loc+'_real" Center="Node" AttributeType="Scalar">\n')
					writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution['fluctuation'][time]['pointData'][key]['magnitude']))+'" Function="$0*(cos($1))">\n')
					writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['magnitude']))+'">\n')
					writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/magnitude\n')
					writer.write('               </DataItem>\n')
					writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['angle']))+'">\n')
					writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/angle\n')
					writer.write('               </DataItem>\n')
					writer.write('           </DataItem>\n')
					writer.write('        </Attribute>\n')
					#Write imaginary part
					writer.write('       <Attribute Name="'+key_loc+'_imaginary" Center="Node" AttributeType="Scalar">\n')
					writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution['fluctuation'][time]['pointData'][key]['magnitude']))+'" Function="$0*sin($1)*(-1)">\n')
					writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['magnitude']))+'">\n')
					writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/magnitude\n')
					writer.write('               </DataItem>\n')
					writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution['fluctuation'][time]['pointData'][key_loc]['angle']))+'">\n')
					writer.write('                   '+solutFileName.rsplit('/')[-1]+':/fluctuation/'+time+'/pointData/'+key_loc+'/angle\n')
					writer.write('               </DataItem>\n')
					writer.write('           </DataItem>\n')
					writer.write('        </Attribute>\n')
				
			writer.write('     </Grid>\n')
			writer.write('  </Domain>\n')
			writer.write('</Xdmf>\n')


