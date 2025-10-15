import  os
import  h5py
import  numpy               as np
from    FELiCS.Fields.Field import Field
from    FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class export:
    """
    This class implements the mapping- and export-Functionality, which is needed
    to visulize the full solution, including the values at the DOF-coordinates.
    Therefore the solution has to be mapped from a P2-Functionspace to a
    P1-Functionspace. Higher dimensional mappings are not possible at the moment.
    This class is a parent of the meanFlowClass and the fluctuationSolutions-class
    in which its methods are used to map and export the meanflow and the
    fluctuations.

    Private Attributes:

    Protected Attributes:
    - _FEMSpaces: Object containing the Cacluation FEM-Spaces, the export
        FEM-Spaces and the corresponding meshes.
    - _param: Object, containing the parameters of the calculation.
    - _exportZeroScalarField: FENiCS scalar Function containing zeros.
    - _exportZeroVectorField: FENics vector-Function containing zeros.
    - _uValuesList: List of strings, containing the names of the velocity
    components.

    """
    def __init__(self, param, FEMSpaces):
        """
        constructor of the export-class

        Function arguments:
        - FEMSpaces: Object containing the Cacluation FEM-Spaces, the export
            FEM-Spaces and the corresponding meshes.
        - param: Object, containing the parameters of the calculation.

        Function returns:

        """
        self._FEMSpaces     = FEMSpaces
        self._param         = param
        self._exportMesh    = FEMSpaces.exportMesh

        self._exportZeroScalarField = Field(self._FEMSpaces.P1Export, self._exportMesh, name="exportZeroScalar")
        self._exportZeroVectorField = Field(
                self._FEMSpaces.FunctionSpaceVectorVelocityExport,
                self._exportMesh, 
                name="exportZeroVector"
                )

    def _mapCalcToExport(self, exportObject):
        """
        this function does the mapping from the calculation space to the export
        Space for the meanflow and the fluctuations. If the given argument is
        of type dictionary, the meanflow-mapping is done. If exportObject is of
        type ndarray, the fluctuation mapping is done.
        The method returns 2 dicts containing Fenics-Functions. exportDictReal
        contains the real part and exportDictImag contains the imaginary part.

        Function arugment:
        - exportObject: meanflow (dict of fenics functions) or fluctuation
        (numpy-array)

        Function Returns:
        - exportDictReal: Dict of Fenics Functions, containing the real part
        - exportDictImag: Dict of Fenics Functions, containing the imaginary part
        """
        from    scipy   import interpolate    
  
        valueDict = {}
        dofsExport = self._exportMesh.coordinates()
        # get the index-vector for the mapping from P2 to P1-Export Space:
        VectorCalcToP1ExportIndecies = self._FEMSpaces.mappingObj.VectorCalcToP1ExportIndecies
        P2CalcToP1ExportIndecies = self._FEMSpaces.mappingObj.P2CalcToP1ExportIndecies

        #VectorCalcToP1VectorExportIndecies = self._FEMSpaces.mappingObj.VectorCalcToP1VectorExportIndecies
        if isinstance(exportObject, dict):
            for indexOfFieldInList, fieldNameFieldToExport in enumerate(
                                                    list(exportObject.keys())
                                                        ):

                ValueArray = exportObject[fieldNameFieldToExport].getCoefficientArray()

                numSubSpaces = exportObject[fieldNameFieldToExport].space.num_sub_spaces
                if numSubSpaces > 1:

                    valueDict[fieldNameFieldToExport] = Field(self._FEMSpaces.FunctionSpaceVectorVelocityExport, self._exportMesh)
                    # dofsCoordsCalc = exportObject[fieldNameFieldToExport].function_space.tabulate_dof_coordinates()
                    #indexVector = self.mappingFunc(dofsCoordsCalc[:, 0:2], dofsExport[:, 0:2])
                    tempSolutionArray = np.zeros((self._param.Case.nDim,dofsExport.shape[0] ), dtype=complex)
                    for subSpaceNum in range(numSubSpaces):
                        indicesOfComponentExport = valueDict[fieldNameFieldToExport].space.sub(subSpaceNum).collapse()[1]
                        indicesOfComponentCalc = self._FEMSpaces.FunctionSpaceVectorVelocity.sub(subSpaceNum).collapse()[1]
                        valueDict[fieldNameFieldToExport].function.x.array[indicesOfComponentExport] = exportObject[fieldNameFieldToExport].function.x.array[indicesOfComponentCalc][P2CalcToP1ExportIndecies]


                else:
                    valueDict[fieldNameFieldToExport] = Field(self._FEMSpaces.P1Export, self._exportMesh, name=fieldNameFieldToExport)
                    #ValueArray = ( flucRealCalc.sub(indexOfFieldInList).collapse().x.array + 1j * flucImagCalc.sub(indexOfFieldInList).collapse().x.array )
                    #subSpaceDofCoordinates = exportObject[fieldNameFieldToExport].function_space.tabulate_dof_coordinates()
                    #indexVector = self.mappingFunc(subSpaceDofCoordinates[:, 0:2], dofsExport[:, 0:2])
                    valueDict[fieldNameFieldToExport].function.x.array[:] = exportObject[fieldNameFieldToExport].function.x.array[P2CalcToP1ExportIndecies]
                    #tempSolutionArray = ValueArray[P2CalcToP1ExportIndecies]
                #valueDict[fieldNameFieldToExport] = tempSolutionArray

                # maybe in the future the interpolation between different meshes
                # will be available again. Then this code can be used.
                # exportDictImag[fieldNameFieldToExport] = zerofield
                # exportDictReal[fieldNameFieldToExport].interpolate(
                #                               exportObject[fieldNameFieldToExport]
                #                                                   )


        elif isinstance(exportObject, np.ndarray):
            # TODO Sophie: give correct names?
            flucRealCalc = Field(self._FEMSpaces.VMixed, None)
            flucImagCalc = Field(self._FEMSpaces.VMixed, None)

            flucRealCalc.function.x.array[:] = np.real(exportObject[:]).astype(
                                                                        float
                                                                            )
            flucImagCalc.function.x.array[:] = np.imag(exportObject[:]).astype(
                                                                        float
                                                                            )

            # linearFunctionReal = Function(self._FEMSpaces.VMixedExport)
            # linearFunctionImag = Function(self._FEMSpaces.VMixedExport)


            dofsExport = self._exportMesh.coordinates()

            for field in self._transportedQuantities:
                if len(self._transportedQuantities) > 1:
                    indexOfFieldInList = self._transportedQuantities\
                    .index(field)
                    numSubSpaces = flucRealCalc.space.sub(indexOfFieldInList).collapse()[0].num_sub_spaces

                    if numSubSpaces > 1:
                        # calculate the complex solution of the vectorfield
                        ValueArray = ( flucRealCalc.function.sub(indexOfFieldInList).collapse().x.array + 1j * flucImagCalc.function.sub(indexOfFieldInList).collapse().x.array )
                        tempSolutionArray = np.zeros((numSubSpaces, self._exportMesh.coordinates().shape[0] ), dtype=complex)
                        # vectorSpaceDofCoords = flucRealCalc.function_space.sub(indexOfFieldInList).collapse()[0].tabulate_dof_coordinates()

                        #indexVector = self.mappingFunc(vectorSpaceDofCoords[:, 0:2], dofsExport[:, 0:2])

                        for subSpaceNum in range(numSubSpaces):

                            # split the complex vectorField in its components
                            # in the ordering of the exportMesh-coordinates

                            # get indices of vectorField-component:
                            dofsOfCalc = flucRealCalc.function.sub(indexOfFieldInList).collapse().function_space.sub(subSpaceNum).collapse()[1]
                            tempSolutionArray[subSpaceNum, :] = ValueArray[dofsOfCalc][VectorCalcToP1ExportIndecies]


                    else:
                        ValueArray = ( flucRealCalc.function.sub(indexOfFieldInList).collapse().x.array + 1j * flucImagCalc.function.sub(indexOfFieldInList).collapse().x.array )
      
                        # Check if the transported quantity was obtained on P2 elts, otherwise we need to 
                        # interpolate from P1 to P2 meshes
                        if len(ValueArray) < len(dofsExport):
                            ##TODO: Sophie: I changed the command, since for the newer dolfinx versions (>0.5.0) the FelicsMesh cannot be given to the FunctionSpace anymore.
                            ## This command should be wrapped in the future.
                            #dofsSolP1 = self._FEMSpaces.P1.mesh.coordinates()
                            meshP1              = self._FEMSpaces.P1.mesh
                            gdim                = meshP1.topology.dim
                            dofsSolP1           = meshP1.geometry.x[:, 0:gdim]
                            tempSolutionArray   = interpolate.griddata(dofsSolP1,ValueArray,dofsExport,method='linear')
                            logger.debug(f'{field}-fluctuations P1-function interpolated on export mesh.')
       
                        else:
                            # subSpaceDofCoordinates = flucRealCalc.sub(indexOfFieldInList).collapse().function_space.tabulate_dof_coordinates()
                            #indexVector = self.mappingFunc(subSpaceDofCoordinates[:, 0:2], dofsExport[:, 0:2])
                            tempSolutionArray = ValueArray[P2CalcToP1ExportIndecies]
                    valueDict[field] = tempSolutionArray
                else:
                    # muss noch angepasst werden!
                    ValueArray = ( flucRealCalc.getCoefficientArray() + 1j * flucImagCalc.getCoefficientArray() )
                    valueDict[field] = ValueArray[P2CalcToP1ExportIndecies]

        else:
            print(f'Cant export the given object of type {type(exportObject)}. Please provide an object of type Dict or np.ndarray!')

        return valueDict


    def _calculateVertexValuesFromDict(
                                        self,
                                        dictReal,
                                        dictImag,
                                        ):
        """
        this function calculates the values of the vertices of the exportMesh
        for each field and brings them in the order of the coordinates of the
        exportMesh.

        function arguments:
        - dictReal: Dictonary of Functions of the the export Space, containing
        the real-values of the field.
        - dictImag: Dictonary of Functions of the the export Space, containing
        the imaginary-values of the field.

        function returns:
        - valueDict: Dictonary containing the vertex-values for each field.

        """
        valueDict = {}
        for key in list(dictReal.keys()):
            tempSolutionArray = ( dictReal[key]\
            .compute_vertex_values() +
                1j * dictImag[key].compute_vertex_values() )
            # The vector components (velocity u) need to be reshaped
            if dictReal[key].function_space().num_sub_spaces() > 1:
                valueDict[key] = tempSolutionArray.reshape(
                    (self._param.Case.getNVelocityComponents(),-1)
                    )
            else:
                valueDict[key] = tempSolutionArray

        return valueDict


    def _createH5GroupStructure(self, hf):
        """
        this function creates the basic group-structure in the opened h5-file.

        function arguments:
        - fileHandlerList:
        -appendFlag: this boolean flag specifies if the fields are appended to
        a given file.

        function returns:
        - meanflowGroupList: list of group-position in h5-file. One element of
        the list points into the meanflow group of a h5-file.
        - pointDataGroupList: list of group-position in h5-file. One element of
        the list points into the pointData group of a h5-file.
        """
        #pointDataGroupList = []
        #for hf in fileHandlerList:
        # create the Group-Structure for fluctuation:
        if 'fluctuation' not in hf.keys():
            fluctuationGroup = hf.create_group('fluctuation')
        else:
            fluctuationGroup = hf['fluctuation']

        # get the number of last exorted solution inside the hf-File:
        lastExportedSolutionList = [int(i)for i in fluctuationGroup.keys()]

        if len(lastExportedSolutionList) == 0:
            currentExportNumberGroup = fluctuationGroup.create_group('0')
        else:
            lastExportedSolutionNumber = max(lastExportedSolutionList)
            currentExportNumberGroup = fluctuationGroup.create_group(f'{lastExportedSolutionNumber+1}')


        if np.imag(np.round(self._omega, 3)) < 0:
            eValStr = f"{np.real(np.round(self._omega, 3))}{np.imag(np.round(self._omega, 3))}j"
        else:
            eValStr = f"{np.real(np.round(self._omega, 3))}+{np.imag(np.round(self._omega, 3))}j"


        if eValStr not in currentExportNumberGroup.keys():
            freqGroup = currentExportNumberGroup.create_group(eValStr)
            pointDataGroup = currentExportNumberGroup.create_group('pointData')

        else:
            freqGroup = currentExportNumberGroup[eValStr]
            pointDataGroup = currentExportNumberGroup['pointData']


        currentExportNumberGroup.attrs.create('analysisType', data=self._param.Case.AnalysisMode)
        currentExportNumberGroup.attrs.create('solutionType', data=self.solutionKind)
        currentExportNumberGroup.attrs.create('frequency', data=eValStr)

        return pointDataGroup

    def _writeDictToH5(
                        self,
                        exportDict,
                        hf,
                        exportAngle = False,
                        ):
        """
        this function exports a given dict to a folder in the H5-File specified
        by hf. If the flag exportAngle is set, then angle and magnitude are
        exported. Else only the magnitude is exported.

        Function Arguments:
        - exportDict: Includes the Fielddata and the Field Names. Either its a
        MeanflowDict or a fluctDict.
        - hf: h5py File-Object. It points on a Folder in the opened H5-File, in
        which the Data should be exported.
        - exportAngle: Optional Flag, which specifies if the angle should be
        exported or  not. The default-value is False.

        Function Returns:
        - hf: Returns the h5py- File Object in which the exportDict was written
        """
        for i, groupName in enumerate(list(exportDict.keys())):
            #if exportDict[groupName].shape[0] == self._param.Case.nDim:
            if groupName in ['u', 'rhou', 'u_forcing_r', 'u_forcing_i']:
                for i, component in enumerate(self._param.BoundaryCondition.VelocityComponents):
                    if 'u_' in groupName:
                        uComponent = 'u' + component
                        newGroupName = groupName.replace('u_', f'{uComponent}_')
                        currFlowVarGroup = hf.create_group(newGroupName)
                    else:
                        uComponent = groupName + component
                        currFlowVarGroup = hf.create_group(uComponent)

                    # if the boolean exportAngle is True, values of typ flucs are
                    # written to h5. Else the meanflow is written. In case of export
                    # of meanflow, no np.abs() is called on the data.
                    if exportAngle:
                        currFlowVarGroup.create_dataset(
                        'magnitude',
                        data = np.abs(exportDict[groupName][i, :]),)
                        currFlowVarGroup.create_dataset(
                        'angle',
                        data = np.angle(exportDict[groupName][i,:]),
                                                )
                    else:
                        currFlowVarGroup.create_dataset(
                        'magnitude',
                        data = np.real(exportDict[groupName][i, :]),)
            else:
                currFlowVarGroup = hf.create_group(groupName)

                # if the boolean exportAngle is True, values of typ flucs are
                # written to h5. Else the meanflow is written. In case of export
                # of meanflow, no np.abs() is called on the data.
                if exportAngle:
                    # print(groupName)
                    # print(type(groupName))
                    # print(type(exportDict[groupName]))
                    # print(type(np.abs(exportDict[groupName])))
                    # print(exportDict[groupName])
                    currFlowVarGroup.create_dataset(
                    'magnitude',
                    data = np.abs(exportDict[groupName]),
                    )
                    currFlowVarGroup.create_dataset(
                    'angle',
                    data = np.angle(exportDict[groupName][:]),
                                            )
                else:
                    currFlowVarGroup.create_dataset(
                    'magnitude',
                    data = np.real(exportDict[groupName]),)
        return hf



    def _XMFwriteHeader(
                        self,
                        writer,
                        time,
                        meshh5,
                        relMeshPath="",
                        ):
        """
        This file writes the begining of the xmf-file containing the Header and
        Mesh information

        Function Input:
        - writer: The Reference to the opened xmf-file
        - time: The corresponding eigenvalue, which should be exported as time
        value into paraview
        - meshh5: File Reference to the opened h5-File, containing the mesh

        Function Returns:

        """

        tri = meshh5['cells']['triangles']
        filename = self._meshfilename.split(self._param.Export.ExportFolder + '/')[1]

        relPathWithFileName = f"{relMeshPath}{filename}"

        if self._param.Case.nDim == 2:
            cellStyle = "Triangle"
        else:
            cellStyle = "Tetrahedron"

        writer.write('<?xml version="1.0" ?>\n')
        writer.write(
    '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n'
                    )
        writer.write('  <Domain>\n')
        writer.write('     <Grid Collection="Triangle_Mesh'+str(time)+'" Name="solution-Triangle'+str(time)+'">\n')
        writer.write('     <Time Value=" '+str(np.real(time))+'" />\n')
        writer.write('        <Topology Type="'+ cellStyle +'" NumberOfElements="     '+str(len(tri[:,0]))+'">\n')
        writer.write('          <DataItem ItemType="Function" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'" Function="$0 - 0">\n')
        writer.write('           <DataItem Format="HDF" DataType="Int" Dimensions="'+str(len(tri[0,:])*len(tri[:,0]))+'">\n')
        writer.write('              '+relPathWithFileName+':/cells/triangles\n')
        writer.write('           </DataItem>\n')
        writer.write('          </DataItem>\n')
        writer.write('        </Topology>\n')
        writer.write('        <Geometry Type="X_Y">\n')
        for key in list(meshh5['coordinates'].keys()):
            writer.write('           <DataItem Format="HDF" ItemType="Uniform" Precision="8" NumberType="Float" Dimensions="     '+str(len(meshh5['coordinates'][key]))+'">\n')
            writer.write('              '+relPathWithFileName+':/coordinates/'+key+'\n')
            writer.write('           </DataItem>\n')
        writer.write('        </Geometry>\n')


    def _XMFwriteMeanflow(
                        self,
                        writer,
                        meanflowFileHandler
                        ):
        """
        this method writes the the part of the xmf-file, corresponding to
        the meanflow.

        Function arguments:
        - writer: The Reference to the opened xmf-file
        - meanflowFileHandler: File Reference to the meanflow .h5-File

        Functions returns:
        -
        """
        for groupName in list(meanflowFileHandler['meanflow'].keys()):
            key = f'meanflow/{groupName}'
            key_loc = groupName
            suffix = 'MeanFlow_'
            #Write magnitude
            writer.write('       <Attribute Name="'+ suffix +key_loc +'" Center="Node" AttributeType="Scalar">\n')
            writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(meanflowFileHandler[key]['magnitude']))+'" Function="$0">\n')
            writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(meanflowFileHandler[key]['magnitude']))+'">\n')
            writer.write('                   '+meanflowFileHandler.filename.rsplit('/')[-1]+':/'+key+'/magnitude\n')
            writer.write('               </DataItem>\n')
            writer.write('           </DataItem>\n')
            writer.write('        </Attribute>\n')



    def _XMFwriteFlucts(
                        self,
                        writer,
                        listOfPaths,
                        listOfKeyloc,
                        solution,
                        solutFileName,
                        ):
        """
        this method writes the part of the xmf-file corresponding to the
        fluctuations

        Function arguments:
        - writer: The Reference to the opened xmf-file
        - listOfPaths: List of the path to the different Fields in the
        fluctuation-.h5-File.
        - listOfKeyloc: List of Strings, representing the Name, which is
        displayed in paraview, i.e. Forcing_ux
        - solution: File Refrence to the opened .h5 file, in which the
        fluctuations of the current fluctsolutObj are stored.
        - solutFileName: String, representing the File Name of the .xmf file.

        Function returns:

        """

        for i, key in enumerate(listOfPaths):
            key_loc = listOfKeyloc[i]
            if 'angle' in list(solution[key].keys()):
                suffix = '_magnitude'
            else:
                suffix = ''
            #Write magnitude
            writer.write('       <Attribute Name="'+key_loc+suffix+'" Center="Node" AttributeType="Scalar">\n')
            writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution[key]['magnitude']))+'" Function="$0">\n')
            writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['magnitude']))+'">\n')
            writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/magnitude\n')
            writer.write('               </DataItem>\n')
            writer.write('           </DataItem>\n')
            writer.write('        </Attribute>\n')
            if 'angle' in list(solution[key].keys()):
                #Write angle
                writer.write('       <Attribute Name="'+key_loc+'_angle" Center="Node" AttributeType="Scalar">\n')
                writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution[key]['angle']))+'" Function="$0">\n')
                writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['angle']))+'">\n')
                writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/angle\n')
                writer.write('               </DataItem>\n')
                writer.write('           </DataItem>\n')
                writer.write('        </Attribute>\n')
                #Write real part
                writer.write('       <Attribute Name="'+key_loc+'_real" Center="Node" AttributeType="Scalar">\n')
                writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution[key]['magnitude']))+'" Function="$0*(cos($1))">\n')
                writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['magnitude']))+'">\n')
                writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/magnitude\n')
                writer.write('               </DataItem>\n')
                writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['angle']))+'">\n')
                writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/angle\n')
                writer.write('               </DataItem>\n')
                writer.write('           </DataItem>\n')
                writer.write('        </Attribute>\n')
                #Write imaginary part
                writer.write('       <Attribute Name="'+key_loc+'_imaginary" Center="Node" AttributeType="Scalar">\n')
                writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution[key]['magnitude']))+'" Function="$0*(sin($1))">\n')
                writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['magnitude']))+'">\n')
                writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/magnitude\n')
                writer.write('               </DataItem>\n')
                writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['angle']))+'">\n')
                writer.write('                   '+solutFileName.rsplit('/')[-1]+':/'+key+'/angle\n')
                writer.write('               </DataItem>\n')
                writer.write('           </DataItem>\n')
                writer.write('        </Attribute>\n')



    def writeXMFFile(
                self,
                meshFile,
                meanflowFilename,
                fieldsFile,
                appendFlag=False,
                    ):
        """
        This method creates a XDMF-File, which holds references between the mesh
        and the data-file. Futhermore with this file a import into paraview is
        simplified.

        Function Arguments:
        - meshFile: String-Name of the h5-file, containing the mesh information
        -fieldsFile: String-Name of the h5-file, containing the field
        information
        - meanflowFilename: String, which represents the Filename, in which the
        meanflow is saved

        Function returns:

        """

        if '/' in fieldsFile.filename:
            outputFolder,outputFile=fieldsFile.filename.rsplit('/', 1)
        else:
            outputFolder='.'
            #outputFile=fieldsFile.filename

        meshh5 = h5py.File(meshFile, 'r')

        self._meshfilename = meshFile

        meanflowFileHandler = h5py.File(f'{self._param.Export.ExportFolder}/{meanflowFilename}', 'r')

        #for i,solution in enumerate(fieldsFile):
        solutFileName = fieldsFile.filename
        solution = fieldsFile

        TimeStamps=list(solution['fluctuation'].keys())

        for i_time,time in enumerate(TimeStamps):
            #pdb.set_trace()
            time = solution[f'fluctuation/{time}'].attrs.get('frequency')
            #time = solution['fluctuation/time'].attrs.get('frequency')
            #writer = open(solutFileName[:-3]+'_'+'.xmf', 'w')
            time = time.replace('+-','+')
            time = time.replace('[','')
            time = time.replace(']','')
            #Sophie: I put this here, in case the string is not a complex number. TODO: handle differently?
            time_str = str(np.real(complex(time)))
            if time_str in solutFileName:
                filename = solutFileName[:-3] + '.xmf'
            else:
                filename = solutFileName[:-3]+'_'+time_str+'.xmf'

            writer = open(filename, 'w')

            # write header and meshinformation to xdmf-file:
            self._XMFwriteHeader(writer, time, meshh5)

            # write Mean Flow Information to Xdmf-File:
            self._XMFwriteMeanflow(writer, meanflowFileHandler)

            listOfPaths = []
            listOfKeyloc = []
            for curSolution in list(solution['fluctuation'].keys()):
                solutionKind = solution[f'fluctuation/{curSolution}'].attrs.get('solutionType')
                for key in list(solution[f'fluctuation/{curSolution}/pointData/'].keys()):
                    listOfPaths.append(f'fluctuation/{curSolution}/pointData/{key}')
                    listOfKeyloc.append(f'{solutionKind}_{key}')

            self._XMFwriteFlucts(writer, listOfPaths, listOfKeyloc, solution, solutFileName)

            writer.write('     </Grid>\n')
            writer.write('  </Domain>\n')
            writer.write('</Xdmf>\n')
            writer.close()

            if self._param.Export.Video:

                # params of the video-export:
                nSnaps = 100
                videoOutputFolderName = 'videoOutput'
                outputFolder, outputFile = solutFileName.rsplit('/', 1)
                basefilename = f'{outputFolder}/{videoOutputFolderName}/{outputFile[:-3]}_video'
                os.chdir(self._param.Export.ExportFolder)
                # create videooutput-Folder, if it not exist.
                if videoOutputFolderName not in os.listdir():
                    os.mkdir(videoOutputFolderName)

                os.chdir("..")
                # angular step size in radiants:
                phiStep = -(360/nSnaps)*np.pi/180

                solutFileName = f'../{outputFile}'
                for snap in range(nSnaps):

                    delPhi = phiStep*snap

                    filename = basefilename + str(snap) + '.xmf'
                    writer = open(f'{filename}', 'w')
                    self._XMFwriteHeader(writer, snap, meshh5, "../")
                    for i, key in enumerate(listOfPaths):
                        key_loc = listOfKeyloc[i]
                        writer.write('       <Attribute Name="'+key_loc+'" Center="Node" AttributeType="Scalar">\n')
                        writer.write('           <DataItem ItemType="Function" Dimensions="'+str(len(solution[key]['magnitude']))+f'" Function="$0*(cos((-1)*($1+({delPhi}))))">\n')
                        writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['magnitude']))+'">\n')
                        writer.write('                   '+solutFileName+':/'+key+'/magnitude\n')
                        writer.write('               </DataItem>\n')
                        writer.write('               <DataItem Format="HDF" Precision="8" NumberType="Float" Dimensions="     '+str(len(solution[key]['angle']))+'">\n')
                        writer.write('                   '+solutFileName+':/'+key+'/angle\n')
                        writer.write('               </DataItem>\n')
                        writer.write('           </DataItem>\n')
                        writer.write('        </Attribute>\n')
                    writer.write('     </Grid>\n')
                    writer.write('  </Domain>\n')
                    writer.write('</Xdmf>\n')
                    writer.close()
