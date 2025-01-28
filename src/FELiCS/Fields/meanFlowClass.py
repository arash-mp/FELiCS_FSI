# Third party libraries
from h5py import File
from dolfinx.fem import Function

# Local libraries and methods
from FELiCS.Misc.functions import printWarning, printDebug
from FELiCS.IO.export import export
from FELiCS.Fields.fieldProperties import fieldProperties
from FELiCS.Equation.dependentVariables.energyHandler import energyHandler
from FELiCS.Equation.dependentVariables.equationOfStateHandler import equationOfStateHandler
from FELiCS.Equation.dependentVariables.reactionHandler import reactionHandler
from FELiCS.Misc.tensorUtils import Tensor

class meanFlowClass(
    fieldProperties,
    energyHandler,
    reactionHandler,
    equationOfStateHandler,
    export,
):
    """
    Parent classes:
    - export
    - fieldProperties

    Child classes:

    """
    def __init__(
            self,
            param,
            FEMSpaces,
            mesh,
    ):

        self._isMean = True
        self._isFluctuation = False
        self._param = param
        self._FEMSpaces = FEMSpaces
        self._mesh = mesh
        self._coordinateSystem = mesh.coordinateSystem
        fieldProperties.__init__(self)
        self._meanflowFilename = None
        self.__mixture = param.Mixture
        self._zeroField = Function(self._FEMSpaces.P2)
        self._zeroFieldTensor = Tensor(
                                       Function(self._FEMSpaces.P2),
                                       self._coordinateSystem,
                                       )
        self._oneFieldArray = Function(self._FEMSpaces.P2)
        self._oneFieldArray.x.array[:] = 1.0
        self._oneField = Tensor(
            self._oneFieldArray,
            self._coordinateSystem, 
            )
        
        #self.addDerivativeFieldsToMean()
        #self.initLamDiff()
        #if param.Case.Reaction:
        #    print('1')
            #print(param.Mixture.getReactionMechanism()['type'])
            #if param.Mixture.getReactionMechanism()['type'] == '2S-SM2':
            #    print('2')
            #    param.Mixture.getReactionMechanism()['type']
            #    self.calculateSpeciesEnthalpy()
            #    from FELiCS.Equation.Reactions.c2sm2 import C2SM2
            #    YCH4_lim = 0.043 * 1e-4
            #    self.__reaction = C2SM2(YCH4_lim, 2)
            #    self.__reaction.computeSensitivities(self.T,
            #                                         self.rho,
            #                                         self.Y('CH4'),
            #                                         self.Y('CO'),
            #                                         self.Y('O2'),
            #                                         self.Y('CO2'))
            #    self._fieldDict['Q'] = self.__reaction.Q(self.T,
            #                                                self.rho,
            #                                                self.Y('CH4'),
            #                                                self.Y('CO'),
            #                                                self.Y('O2'),
            #                                                self.Y('CO2'))
            #if param.Mixture.ReactionMechanism['type'] == 'NOx':
            #    print('3')
            #    from FELiCS.Equation.Reactions.NOx import NOx
            #    self.__reaction = NOx(2)

    def importDataFromFile(self):
        printDebug(True,'--------------------------------')
        printDebug(True,'-- Reading InputFlow...')

        self._fieldDict = {}
        self.__notInFileList = []
        self.__RawFlowDict = {}
        self.__nDimRawData = 0
        self._ScalarFunctionSpace = self._FEMSpaces.P2
        self.__VectorFunctionSpace = self._FEMSpaces.FunctionSpaceVectorVelocity
        self.__CoordinateSystemInputData = 'Unknown'
        # Check which type the input file is and read
        if self._param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'h5':
            # In a hdf5 file the data is already interpolated on the mesh
            # from FELiCS.Import import importHDF5File
            self.importHDF5File2()
        else:
            if self._param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'fel':
                self.importFelicsFile()
            if self._param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'xdmf':
                self.importXDMFFile()
            if self._param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'mat':
                self.importMatFile()
            # If necessary perform coordinate transformation (so far only from
            # cartesian to cylindrical)
            self.CoordinateTransformation()
            self.InterpolateOnFELiCSMesh()
            del self.__RawFlowDict
            self.raiseNotInFileListWarning()
        self.initLamDiff()

#    def addDerivativeFieldsToMean(self):
#        from ufl import sqrt
#        from dolfinx.fem import Function
#
#        if self._param.Case.MolViscModel == 'Constant':
#            self._fieldDict['nulam'], \
#            self._fieldDict['nulam'].x.array[:] = self.getConstVisc()
#        #if self._param.Case.Compressible:
#        #    self._fieldDict['c'].interpolate(sqrt(self.gamma * self.p / self.rho))
#
    def initLamDiff(self):
        from dolfinx.fem import Function

        
        if self._param.Case.MolViscModel == 'Constant':
            self._fieldDict['nulam'] = Function(self._FEMSpaces.P2)
            self._fieldDict['nulam'].x.array[:] = self._param.Case.MolVisc

        for specie in self._param.Mixture.getSpeciesList('transported'):
            Sc = self._param.Mixture.species[specie]['Sc']
            nuTot = Function(self._ScalarFunctionSpace)

            if 'nulam' in list(self._fieldDict.keys()):
                nuTot.x.array[:] += self._fieldDict['nulam'].x.array[:]
            if 'nuturb' in list(self._fieldDict.keys()):
                nuTot.x.array[:] += self._fieldDict['nuturb'].x.array[:]
            if 'nuSGS' in list(self._fieldDict.keys()):
                nuTot.x.array[:] += self._fieldDict['nuSGS'].x.array[:]
            self._fieldDict['D_' + specie] \
                = Function(self._FEMSpaces.P2)
            self._fieldDict['D_' + specie].vector[:] \
                = nuTot.x.array[:] / Sc

    def importMatFile(self):
        import scipy.io as spio
        from FELiCS.Misc.functions import printDebug
        filePath = self._param.FlowInput.MeanFlowFilePath
        mat = spio.loadmat(filePath[:-4])
        printDebug(self._param.debug, 'The fields in the Matlab file are '
                   + str(mat.keys()))
        for key in list(mat.keys()):
            self.__RawFlowDict[key] = mat[key]
            self.__RawFlowDict = mat

    def mapToExportMeshAndExport(self, FEMSpaces, filename):
        self._meanflowFilename = filename
        filehandler = File(f'{self._param.Export.ExportFolder}/{filename}', 'w')
        group = filehandler.create_group('meanflow')
        export.__init__(self, self._param, self._FEMSpaces)
        #self._meanfieldDict, dictImag = self._mapCalcToExport(self._fieldDict)
        self._meanfieldDict = self._mapCalcToExport(self._fieldDict)
        #exportDict = self._calculateVertexValuesFromDict(self._meanfieldDict,
        #                                                 dictImag)
        #self._writeDictToH5(exportDict, group)
        self._writeDictToH5(
            self.getVertexValues()._fieldDict,
            group,
        )

    def exportBaseFlowAsHDF5(self, meanflowFilename = 'meanflow.h5'):

        #from fenics import HDF5File
        #mesh = self.fieldDict[list(self.fieldDict.keys())[0]].\
        #    function_space().mesh()
        #exportFilePath = self._param.FlowInput.MeanFlowFilePath[0:-4] + "_" \
        #                 + self._param.Case.MeshFilePath[0:-3].split('/')[-1] \
        #                 + "hdf5"
        # export the mapped Meanflow:
        #meanflowFilename = 'meanflow.h5'
        self.mapToExportMeshAndExport(self._FEMSpaces, meanflowFilename)
        #hdf5file = HDF5File(mesh.mpi_comm(), exportFilePath, 'w')
        #for name in self.fieldDict.keys():
        #    hdf5file.write(self.fieldDict[name], name)

    def importHDF5File2(self):
        from dolfinx.fem import Function
        import h5py
        import numpy as np
        meanflowH5 = h5py.File(f"meanflow.h5", 'r')
        exportMeshH5 = h5py.File(f"{self._param.Case.AnalysisMode}_mesh.h5", 'r')
        coordNameList = ['x', 'y', 'z']

        coordArray = np.zeros(
            (exportMeshH5['coordinates/x'][:].shape[0], self._param.Case.nDim))
        velocityArray = np.zeros(
            (exportMeshH5['coordinates/x'][:].shape[0], self._param.Case.nDim))

        coordArray[:, 0] = exportMeshH5['coordinates/x'][:]
        for i in range(self._param.Case.nDim - 1):
            coordArray[:, i + 1] = exportMeshH5[
                                       f'coordinates/{coordNameList[i + 1]}'][:]

        velocityComponents = self._param.getInternalVelocityComponents()
        for index, comp in enumerate(velocityComponents):
            velocityComponents[index] = f'u{comp}'

        def mappingFunc(exportMeshDOFCoordinates, calcMeshDOFCoordinates):
            # mapping from CalcMesh to exportMesh needs to be done:
            mapping = np.zeros(exportMeshDOFCoordinates.shape[0], dtype=int)
            for index, coordinate in enumerate(calcMeshDOFCoordinates):
                mapping[index] = np.where(
                    np.isclose(coordinate, exportMeshDOFCoordinates).all(
                        axis=1) == True)[0]
            return mapping

        coordinatesOfP2Mesh = self._FEMSpaces.P2.tabulate_dof_coordinates()[:,
                              0:self._param.Case.nDim]

        indexMappingArray = mappingFunc(coordArray, coordinatesOfP2Mesh)

        fieldDict = {}
        nameListMean = self._getMeanFieldsToBeRead()    
        for name in nameListMean:
            if name[0] == 'u' and not (name == 'ut' or name == 'ut_forcing'):
                fieldDict[name] = Function(
                    self._FEMSpaces.FunctionSpaceVectorVelocity)
            else:
                fieldDict[name] = Function(self._FEMSpaces.P2)

            if name == 'u':
                for index, component in enumerate(velocityComponents):
                    indicesOfSubField = \
                    fieldDict[name].function_space.sub(index).collapse()[1]
                    fieldDict[name].x.array[indicesOfSubField] = \
                    meanflowH5[f'meanflow/{component}/magnitude'][:][
                        indexMappingArray]

            else:
                if name in list(meanflowH5[f'meanflow'].keys()):
                    fieldDict[name].vector[:] = \
                    meanflowH5[f'meanflow/{component}/magnitude'][:][
                        indexMappingArray]
                else:
                    self.__notInFileList.append(name)
        if 'ut' in list(fieldDict.keys()):
            fieldDict['ut'].x.array[:] = 0.0
        self._fieldDict = fieldDict
        #for key in list(fieldDict.keys()):
        #    self._fieldDict[key] = Tensor(
        #                            fieldDict[key],
        #                            self.__coordSys,
        #                            )

    def importFelicsFile(self):
        import h5py
        import numpy as np
        Case = self._param.Case
        FlowInput = self._param.FlowInput
        # RawFlowDict is the dictionary directly loaded from the input file
        self.__RawFlowDict = {}
        # List of fields, which are not in the import file
        self.__notInFileList = []
        ## Get Mean flow names
        nameListMean = self._getMeanFieldsToBeRead()
        ## Get AVBP mesh file path
        filePath = self._param.FlowInput.MeanFlowFilePath
        # Open hdf5 file
        h5file = h5py.File(filePath, 'r')

        # Analyze dimension and coordinate system of input data based on the
        # available coordinates. Copy the respective coordinates to the
        # self.__RawFlowDict at the same time
        self.__nDimRawData = 0
        self.__CoordinateSystemInputData = 'Unknown'
        if 'x' in list(h5file['MeanFlow'].keys()):
            self.__RawFlowDict['x'] = np.array(h5file['MeanFlow']['x'])
            self.__nDimRawData += 1
        if 'y' in list(h5file['MeanFlow'].keys()):
            self.__RawFlowDict['y'] = np.array(h5file['MeanFlow']['y'])
            self.__nDimRawData += 1
            self.__CoordinateSystemInputData = 'Cartesian'
        if 'z' in list(h5file['MeanFlow'].keys()):
            self.__RawFlowDict['z'] = np.array(h5file['MeanFlow']['z'])
            self.__nDimRawData += 1
            self.__CoordinateSystemInputData = 'Cartesian'
        if 't' in list(h5file['MeanFlow'].keys()):
            self.__RawFlowDict['t'] = np.array(h5file['MeanFlow']['t'])

            self.__nDimRawData += 1
            if self.__CoordinateSystemInputData == 'Cartesian':
                printWarning('Ambiguous input data: both cylindrical \
                             coordinates and cartesian coordinates \
                             present. Check input data!')
            else:
                self.__CoordinateSystemInputData = 'Cylindrical'
        if 'r' in list(h5file['MeanFlow'].keys()):
            self.__RawFlowDict['r'] = np.array(h5file['MeanFlow']['r'])
            self.__nDimRawData += 1
            if self.__CoordinateSystemInputData == 'Cartesian':
                printWarning('Ambiguous input data: both cylindrical \
                             coordinates and cartesian coordinates \
                             present. Check input data!')
            else:
                self.__CoordinateSystemInputData = 'Cylindrical'

        # Copy all remaining fields to the self.__RawFlowDict
        for name in nameListMean:
            if name[0] == 'u':
                for Component in self._param.BoundaryCondition.VelocityComponents:
                    nameComponent = name[:1] + Component + name[1:]
                    if nameComponent in list(h5file['MeanFlow'].keys()):
                        self.__RawFlowDict[nameComponent] \
                            = np.array(h5file['MeanFlow'][nameComponent])
                    else:
                        self.__notInFileList.append(nameComponent)
            else:
                if name in list(h5file['MeanFlow'].keys()):
                    self.__RawFlowDict[name] \
                        = np.array(h5file['MeanFlow'][name])
                else:
                    self.__notInFileList.append(name)
        # Check if dimensions of RawFlowDict are OK, if not, correct it
        for key in list(self.__RawFlowDict.keys()):
            if len(np.shape(self.__RawFlowDict[key])) > 1:
                self.__RawFlowDict[key] = np.squeeze(self.__RawFlowDict[key])
                
    def importXDMFFile(self):
        '''
        Can't give information about CoordinateSystemInputData.
        Instead it will read the coordinate system in case setting.
        Only tested reading 2-D fields.
        '''
        import h5py
        import numpy as np
        Case = self._param.Case
        FlowInput = self._param.FlowInput
        # RawFlowDict is the dictionary directly loaded from the input file
        self.__RawFlowDict = {}
        # List of fields, which are not in the import file
        self.__notInFileList = []
        ## Get Mean flow names
        nameListMean = self._getMeanFieldsToBeRead()
        ## Get AVBP mesh file path
        filePath = self._param.FlowInput.MeanFlowFilePath[:-4] + 'h5'
        # Open hdf5 file
        def hdf5_to_dict(group):
            result = {}
            for key, item in group.items():
                if isinstance(item, h5py.Dataset):
                    result[key] = item[()]  # Add dataset's value to the result dictionary
                elif isinstance(item, h5py.Group):
                    result[key] = hdf5_to_dict(item)  # Recursively call the function for subgroups
            return result

        with h5py.File(filePath, 'r') as h5file:
            Coordinate = h5file['Mesh']['mesh']['geometry'][:]
            field = hdf5_to_dict(h5file['Function'])

        # Analyze dimension and coordinate system of input data based on the
        # available coordinates. Copy the respective coordinates to the
        # self.__RawFlowDict at the same time
        self.__nDimRawData = np.shape(Coordinate)[1]
        if np.shape(Coordinate)[1]>=1:
            self.__RawFlowDict['x'] = Coordinate[:,0]
        if np.shape(Coordinate)[1]>=2:
            if self._param.Case.CoordinateSystem in ['Cartesian']:
                self.__RawFlowDict['y'] = Coordinate[:,1]
            elif self._param.Case.CoordinateSystem in ['Cylindrical']:
                self.__RawFlowDict['r'] = Coordinate[:,1]
        if np.shape(Coordinate)[1]>=3:
            self.__RawFlowDict['z'] = Coordinate[:,1]

        
        
        # Copy all remaining fields to the self.__RawFlowDict
        for name in nameListMean:
            if name[0] == 'u':
                count = 0
                for Component in self._param.BoundaryCondition.VelocityComponents:
                    nameComponent = name[:1] + Component + name[1:]
                    if 'u' in list(field.keys()):
                        self.__RawFlowDict[nameComponent] \
                            = np.array(field['u']['0'][:,count])
                    elif 'real_u' in list(field.keys()):
                        self.__RawFlowDict[nameComponent] \
                            = np.array(field['real_u']['0'][:,count])
                    else:
                        self.__notInFileList.append(nameComponent)
                    count += 1
                del count
            else:
                if name in list(field.keys()):
                    self.__RawFlowDict[name] \
                        = np.array(field[name]['0'])
                elif 'real_' + name in field.keys():
                    self.__RawFlowDict[name] \
                        = np.array(field['real_'+name]['0'])
                else:
                    self.__notInFileList.append(name)
        # Check if dimensions of RawFlowDict are OK, if not, correct it
        for key in list(self.__RawFlowDict.keys()):
            if len(np.shape(self.__RawFlowDict[key])) > 1:
                self.__RawFlowDict[key] = np.squeeze(self.__RawFlowDict[key])

    def InterpolateOnFELiCSMesh(self):

        import numpy as np
        from dolfinx.fem import Function
        from scipy import interpolate
        from FELiCS.IO.Import import ExpandForAverage, ContractAfterAverage

        printDebug(True, '-- Interpolating on FELiCS mesh...')
        nameListMean = self._getMeanFieldsToBeRead()
        # Get mesh data
        mesh = self._ScalarFunctionSpace.mesh
        # For three-dimensional databases restrict domain to reduce the number
        # of basis points and accelerate the interpolation
        if self.__nDimRawData > 2:
            Bound = True
        # For 2D flows this most often is not necessary
        else:
            Bound = False
        Bound = False
        if Bound:
            Xcoords = [k[0] for k in mesh.coordinates()]
            Ycoords = [k[1] for k in mesh.coordinates()]
            boundingBox = [min(Xcoords), max(Xcoords), min(Ycoords),
                           max(Ycoords), -mesh.hmax(), mesh.hmax()]

            X = self.__RawFlowDict['x']
            R = self.__RawFlowDict['r']

            Xidx = np.argwhere((X > boundingBox[0]) &
                               (X < boundingBox[1])).flatten()
            Ridx = np.argwhere((R < boundingBox[3])).flatten()

            IDX = np.intersect1d(Xidx, Ridx)
            x1_mean = np.array(self.__RawFlowDict['x'][IDX])
            x2_mean = np.array(self.__RawFlowDict['y'][IDX])
            x3_mean = np.array(self.__RawFlowDict['z'][IDX])
            points = np.vstack((x1_mean, x2_mean, x3_mean)).T
            for k in self.__RawFlowDict.keys():
                self.__RawFlowDict[k] = self.__RawFlowDict[k][IDX]
        else:
            x1_mean = np.array(self.__RawFlowDict['x'])
            x1_mean = x1_mean.tolist()
            if self._param.Case.CoordinateSystem in ['Cartesian']:
                x2_mean = np.array(self.__RawFlowDict['y'])
            elif self._param.Case.CoordinateSystem in ['Cylindrical']:
                x2_mean = np.array(self.__RawFlowDict['r'])
            x2_mean = x2_mean.tolist()
            if self.__nDimRawData > 2:
                x3_mean = np.array(self.__RawFlowDict['z'])
                x3_mean = x3_mean.tolist()
                points = np.vstack((x1_mean, x2_mean, x3_mean)).T
            else:
                points = np.vstack((x1_mean, x2_mean)).T

        # get coordinates of FELiCS mesh
        dof_coordinatesP2 = self._ScalarFunctionSpace.tabulate_dof_coordinates()
        # In case of a cylindrical coordinate system in the FELiCS calculation
        # and a 3D input flow, azimuthal averaging must be performed
        if self._param.FlowInput.AveragingDirection == 'Azimuthal':
            dof_InterpolationFELiCSMeshP2 = ExpandForAverage(dof_coordinatesP2,
                                                             self._param)
        # For other cases the interpolation points are identical to the FELiCS
        # mesh coordinates
        else:
            dof_InterpolationFELiCSMeshP2 = dof_coordinatesP2
        namesP2 = []

        # Get number of fields to interpolate
        nFieldsToInterpolate = len(nameListMean) \
                               + self._param.BoundaryCondition.nVelocityComponents - 1
        if self._param.Case.AnalysisMode in ['Input-Output']:
            nFieldsToInterpolate += 2 * (self._param.BoundaryCondition.nVelocityComponents - 1)
        nFieldsToInterpolate -= len(self.__notInFileList)
        # Define vmatrix for interpolation basis values, valsP2
        valsP2 = np.zeros((len(
            self.__RawFlowDict[list(self.__RawFlowDict.keys())[0]]),
                           nFieldsToInterpolate + 6))

        m = 0

        # Fill valsP2 with raw data
        for name in nameListMean:
            if name[0] == 'u' and name not in ['ut_forcing_r',
                                               'ut_forcing_i']:
                self._fieldDict[name] = Function(self.__VectorFunctionSpace)
                # All inplane velocity components are defined as vectors.
                # Therefore, for these, iterate through the components
                for component in self._param.BoundaryCondition.VelocityComponents:
                    # Get the name in plus component
                    nameComponent = name[:1] + component + name[1:]
                    if nameComponent in list(self.__RawFlowDict.keys()):
                        valsP2[:, m] \
                            = np.array(self.__RawFlowDict[nameComponent])
                        namesP2.append(nameComponent)
                        m += 1
            # Do the same as above, for all scalars. Here no iteration through
            # components is necessary
            else:
                self._fieldDict[name] = Function(self._ScalarFunctionSpace)
                if name in list(self.__RawFlowDict.keys()):
                    valsP2[:, m] = np.array(self.__RawFlowDict[name])
                    namesP2.append(name)
                    # increment m
                    m += 1

        # Perform interpolation nearest (imprecise)
        temp_vecP2nearest = interpolate.griddata(points, valsP2,dof_InterpolationFELiCSMeshP2[:, 0:self._param.Case.nDim], method='nearest')
        # Perform interpolation nearest (more precise)
        try:
            temp_vecP2 = interpolate.griddata(points, valsP2,
                                              dof_InterpolationFELiCSMeshP2[:, 0:self._param.Case.nDim],
                                              method='linear')
            temp_vecP2[np.isnan(temp_vecP2)] \
                = temp_vecP2nearest[np.isnan(temp_vecP2)]
        except:
            printWarning("Linear interpolation failed... Continuing with \
                         nearest interpolation. This may cause strong \
                         inaccuracies!")
            temp_vecP2 = np.array(temp_vecP2nearest)
        # In case of nan values in the linear interpolation result, fill the nan
        # values with the results from the nearest interpolation
        temp_vecP2[np.isnan(temp_vecP2)] \
            = temp_vecP2nearest[np.isnan(temp_vecP2)]
        # Now write the interpolated Values to the self._fieldDict
        m = 0
        # Do that for every entry in nameListMean
        for name in nameListMean:
            if name[0] == 'u' and name not in ['ut_forcing_r',
                                               'ut_forcing_i']:
                # All inplane velocity components are defined as vectors.
                # Therefore, for these, iterate through the components
                for idx, component in enumerate(self._param.BoundaryCondition.VelocityComponents):
                    nameComponent = name[:1] + component + name[1:]
                    # if not nameComponent in self.__notInFileList:
                    if nameComponent in list(self.__RawFlowDict.keys()):
                        # Get the indices of the components entries in the
                        # vector
                        dofIDX = self.__VectorFunctionSpace.sub(idx).collapse()[1]
                        if self._param.FlowInput.AveragingDirection \
                                == 'Azimuthal':
                            # For this case a azimuthal average is performed by
                            # using the function ContractFromAzimuthalAverage
                            self._fieldDict[name].sub(idx).x.array[dofIDX] \
                                = ContractAfterAverage(dof_coordinatesP2,
                                                       self._param,
                                                       temp_vecP2[:, m])
                        else:
                            # In this case a simple copy of the interpolation
                            # results is sufficient
                            self._fieldDict[name].sub(idx).x.array[dofIDX] \
                                = np.array(temp_vecP2[:, m])
                        # Increment m
                        m += 1
            else:
                if (not name in self.__notInFileList) or (
                        name in ['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt',
                                 'rstrt']):
                    # Do the same as above, for all scalars. Here no iteration
                    # through components is necessary
                    if name in list(self.__RawFlowDict.keys()):
                        if self._param.FlowInput.AveragingDirection \
                                == 'Azimuthal':
                            self._fieldDict[name].vector[:] \
                                = ContractAfterAverage(dof_coordinatesP2,
                                                       self._param,
                                                       temp_vecP2[:, m])
                        else:
                            self._fieldDict[name].vector[:] \
                                = np.array(temp_vecP2[:, m])
                        m += 1

    def raiseNotInFileListWarning(self):
        from FELiCS.Misc.functions import printWarning
        for name in self.__notInFileList:
            printWarning("  -- Field " + name
                         + " not found in the import file! Assuming Field is \
                           zero...")

    def CoordinateTransformation(self):
        import numpy as np
        """ 
        This function checks if a coordinate transform is necessary to make the 
        input data compatible with the FELiCS mesh. 
        If so the transformation will be performed.
        For all cases first the basis coordinate system is shifted to the axis 
        of the Felics mesh (which is by default the x-axis).
        For the moment only Cart2Cyl is implemented and only calculates the 
        cylindrical velocity components from the carthesian ones.
        
        First calculate the rotation matrix R based on the two vectors A and B.
        A (the vector of origin) is given by the user and corresponds to the 
        axis of rotational symmetry or to the symmetry plane in 2D carthesian 
        meshes, which must intersect the origin, currently B is the vector of 
        destination and is currently by default the x axis (The axis or plane 
        of symmetry for the linear analysis is always the x axis)
        
        Edit:
        User is now required to set the axis mentioned above within the gui 
        settings. A dropdown menu is provided were the user can choose from x,
        y, z for both,  the destination(B) and the origin axis(A)
        """

        if (self.__nDimRawData == 3
                and self._param.Case.CoordinateSystem == 'Cylindrical'
                and self.__CoordinateSystemInputData == 'Cartesian'):
            # if origin axis is y or z the rotation matrix R is formed
            # according to choice
            A = np.zeros(3)
            if originAxis == 'x':
                A[0] = 1
            elif originAxis == 'y':
                A[1] = 1
            elif originAxis == 'z':
                A[2] = 1
            B = np.zeros(3)
            B[0] = 1

            if np.all(A - B == 0):
                print("Origin vector and target vector are identical. Skipping \
                      rotation...")
                [self.__RawFlowDict['x'],
                 self.__RawFlowDict['y'],
                 self.__RawFlowDict['z']] \
                    = [self.__RawFlowDict['x'],
                       self.__RawFlowDict['y'],
                       self.__RawFlowDict['z']]
            else:
                print("Rotating the coordinate system to align with the Felics \
                      mesh...")
                v = np.cross(A, B)
                s = np.linalg.norm(v)
                c = np.dot(A, B)
                v_matrix = np.zeros((3, 3))
                v_matrix[1, 0] = v[2]
                v_matrix[2, 0] = -v[1]
                v_matrix[2, 1] = v[0]
                v_matrix[0, 1] = -v[2]
                v_matrix[0, 2] = v[1]
                v_matrix[1, 2] = -v[0]
                print("Rotating base flow mesh to felics mesh...")
                print("Rotation from (" + str(A[0]) + "," + str(A[1]) + ","
                      + str(A[2]) + ") to (" + str(B[0]) + "," + str(B[1]) + ","
                      + str(B[2]) + ")")
                R = np.identity(3) + v_matrix + np.dot(v_matrix, v_matrix) \
                    * (1 - c) / s ** 2
                # Multiply R on the raw coordinates (x,y,z) to obtain the final
                # coordinates
                [self.__RawFlowDict['x'],
                 self.__RawFlowDict['y'],
                 self.__RawFlowDict['z']] \
                    = list(np.dot(R, np.array([self.__RawFlowDict['x'],
                                               self.__RawFlowDict['y'],
                                               self.__RawFlowDict['z']])))
            print("Performing Coordinate Transform from cartesian to \
                  cylindrical coordinates...")
            # Rotate velocities (ux,uy,uz) to the new coordinate system via
            # rotation matrix R, so they become
            if np.all(A - B == 0):
                [self.__RawFlowDict['ux'],
                 self.__RawFlowDict['uy'],
                 self.__RawFlowDict['uz']] \
                    = [self.__RawFlowDict['ux'],
                       self.__RawFlowDict['uy'],
                       self.__RawFlowDict['uz']]
            else:
                [self.__RawFlowDict['ux'],
                 self.__RawFlowDict['uy'],
                 self.__RawFlowDict['uz']] \
                    = list(np.dot(R, np.array([self.__RawFlowDict['ux'],
                                               self.__RawFlowDict['uy'],
                                               self.__RawFlowDict['uz']])))

            if transformType == 'Cart2Cyl':

                # Calculate the radial and azimuthal coordinate
                self.__RawFlowDict['r'] \
                    = (self.__RawFlowDict['z'] ** 2
                       + self.__RawFlowDict['y'] ** 2) ** 0.5
                self.__RawFlowDict['theta'] \
                    = np.arctan2(self.__RawFlowDict['z'],
                                 self.__RawFlowDict['y'])
                # Obtain radial and tangential velocity (cylindrical
                # coordinates) from the angle theta and the velocities ux and
                # uy (cartesian coordinates)
                self.__RawFlowDict['ur'] \
                    = np.cos(self.__RawFlowDict['theta']) \
                      * self.__RawFlowDict['uy'] \
                      + np.sin(self.__RawFlowDict['theta']) \
                      * self.__RawFlowDict['uz']
                self.__RawFlowDict['ut'] \
                    = -np.sin(self.__RawFlowDict['theta']) \
                      * self.__RawFlowDict['uy'] \
                      + np.cos(self.__RawFlowDict['theta']) \
                      * self.__RawFlowDict['uz']

                if 'rstxx' in list(self.__RawFlowDict.keys()):
                    self.__RawFlowDict['rstrr'] \
                        = np.cos(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstyy'] \
                          + np.sin(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstzz'] \
                          + 2 * np.sin(self.__RawFlowDict['theta']) \
                          * np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstyz']
                    self.__RawFlowDict['rsttt'] \
                        = np.cos(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstzz'] \
                          + np.sin(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstyy'] \
                          - 2 * np.sin(self.__RawFlowDict['theta']) \
                          * np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstyz']
                    self.__RawFlowDict['rstxr'] \
                        = np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstxy'] \
                          + np.sin(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstxz']
                    self.__RawFlowDict['rstxt'] \
                        = -np.sin(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstxy'] \
                          + np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstxz']
                    self.__RawFlowDict['rstrt'] \
                        = np.cos(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstyz'] \
                          - np.sin(self.__RawFlowDict['theta']) ** 2 \
                          * self.__RawFlowDict['rstyz'] \
                          + np.sin(self.__RawFlowDict['theta']) \
                          * np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstzz'] \
                          - np.sin(self.__RawFlowDict['theta']) \
                          * np.cos(self.__RawFlowDict['theta']) \
                          * self.__RawFlowDict['rstyy']

                    del self.__RawFlowDict['rstyy'], \
                        self.__RawFlowDict['rstzz'], \
                        self.__RawFlowDict['rstxy'], \
                        self.__RawFlowDict['rstxz'], \
                        self.__RawFlowDict['rstyz']

    def plot(self, field='all'):
        import matplotlib.pyplot as plt
        from fenics import plot
        import matplotlib as mpl
        import numpy as np
        mpl.use('TkAgg')
        MeanFieldList = list(self._fieldDict.keys())
        if field == 'all':
            PlotList = MeanFieldList
            PlotListTemp = []
            for entry in PlotList:
                if entry[0] == 'u':
                    index = PlotList.index(entry)
                    PlotListTemp.append('ux' + entry[1:])
                    PlotListTemp.append('uy' + entry[1:])
                else:
                    PlotListTemp.append(entry)
            PlotList = PlotListTemp
        else:
            if field in MeanFieldList or field[0] == 'u':
                PlotList = [field]
            else:
                raise Exception(field + ' is not listed in the mean field...')
        for pfield in PlotList:
            if pfield[0] == 'u':
                if pfield[1] == 'x':
                    index = 0
                elif pfield[1] in ['y','r']:
                    index = 1
                elif pfield[1] == 't':
                    index = 2
                fieldToPlot = self._fieldDict[pfield[0] + pfield[2:]][index]
                dofIDX = self._fieldDict[pfield[0] + pfield[2:]].\
                    function_space().sub(index).dofmap().dofs()
                cbarmin = np.min(self._fieldDict[pfield[0] + pfield[2:]].
                                 split()[index].vector()[dofIDX])
                cbarmax = np.max(self._fieldDict[pfield[0] + pfield[2:]].
                                 split()[index].vector()[dofIDX])
            else:
                fieldToPlot = self._fieldDict[pfield]
                cbarmin = np.min(self._fieldDict[pfield].vector()[:])
                cbarmax = np.max(self._fieldDict[pfield].vector()[:])

            fig = plt.figure()
            ax = fig.add_subplot(2, 1, 1)
            cs = plot(fieldToPlot)
            plt.title(pfield)
            ax = fig.add_subplot(2, 1, 2)
            norm = mpl.colors.Normalize(vmin=cbarmin, vmax=cbarmax)
            cb1 = mpl.colorbar.ColorbarBase(ax,
                                            ticks=[cbarmin, cbarmax],
                                            norm=norm,
                                            orientation='horizontal')
            # cbar = fig.colorbar(cs,ticks=[cbarmin,cbarmax])
            # cbar.set_clim(-cbarmin, cbarmax)
            plt.show()

    def calculateSpeciesEnthalpy(self):
        from fenics import project
        # noinspection PyUnresolvedReferences
        import FELiCS.Equation.Reactions.janafspecie as janafspecie
        # noinspection PyUnresolvedReferences
        import FELiCS.Equation.Reactions.janafopenfoam as janafopenfoam
        # initialize janafOF
        janaf = janafopenfoam.Janafopenfoam(2)

        N2 = janafspecie.Janafspecie(28.0134, 200, 5000, 1000,
                                     [3.29868, 0.00140824, -3.96322e-06,
                                      5.64152e-09, -2.44486e-12, -1020.9,
                                      3.95037],
                                     [2.92664, 0.00148798, -5.68476e-07,
                                      1.0097e-10, -6.75335e-15, -922.798,
                                      5.98053],
                                     2)

        CO2 = janafspecie.Janafspecie(44.01, 200, 5000, 1000,
                                      [2.27572, 0.00992207, -1.04091e-05,
                                       6.86669e-09, -2.11728e-12, -48373.1,
                                       10.1885],
                                      [4.45362, 0.00314017, -1.27841e-06,
                                       2.394e-10, -1.66903e-14, -48967,
                                       -0.955396],
                                      2)

        O2 = janafspecie.Janafspecie(31.9988, 200, 5000, 1000,
                                     [3.21294, 0.00112749, -5.75615e-07,
                                      1.31388e-09, -8.76855e-13, -1005.25,
                                      6.03474],
                                     [3.69758, 0.00061352, -1.25884e-07,
                                      1.77528e-11, -1.13644e-15, -1233.93,
                                      3.18917],
                                     2)

        CH4 = janafspecie.Janafspecie(16.043, 200, 5000, 1000,
                                      [0.778741, 0.0174767, -2.78341e-05,
                                       3.04971e-08, -1.22393e-11, -9825.23,
                                       13.7222],
                                      [1.68348, 0.0102372, -3.87513e-06,
                                       6.78559e-10, -4.50342e-14, -10080.8,
                                       9.6234],
                                      2)

        CO = janafspecie.Janafspecie(28.0106, 200, 5000, 1000,
                                     [3.26245, 0.00151194, -3.88176e-06,
                                      5.58194e-09, -2.47495e-12, -14310.5,
                                      4.8489],
                                     [3.02508, 0.00144269, -5.63083e-07,
                                      1.01858e-10, -6.91095e-15, -14268.4,
                                      6.10822],
                                     2)
        # H2O gaseous
        H2O = janafspecie.Janafspecie(18.0153, 200, 5000, 1000,
                                      [3.38684, 0.00347498, -6.3547e-06,
                                       6.96858e-09, -2.50659e-12, -30208.1,
                                       2.59023],
                                      [2.67215, 0.00305629, -8.73026e-07,
                                       1.201e-10, -6.39162e-15, -29899.2,
                                       6.86282],
                                      2)

        self.__hSpec = {}
        self.__hSpec['CH4'] = project(janaf.janaf_Hs_expr(CH4, self.T),
                                      self._FEMSpaces.P2)
        self.__hSpec['O2'] = project(janaf.janaf_Hs_expr(O2, self.T),
                                     self._FEMSpaces.P2)
        self.__hSpec['CO'] = project(janaf.janaf_Hs_expr(CO, self.T),
                                     self._FEMSpaces.P2)
        self.__hSpec['CO2'] = project(janaf.janaf_Hs_expr(CO2, self.T),
                                      self._FEMSpaces.P2)
        self.__hSpec['H2O'] = project(janaf.janaf_Hs_expr(H2O, self.T),
                                      self._FEMSpaces.P2)

    def smoothField(self, field_name, n_timesteps, dt, alpha):
        from fenics import (
            Constant,
            dot,
            grad,
            dx,
            TrialFunction,
            TestFunction,
            Function,
            lhs,
            rhs,
            solve,
            File
        )
        field = self._fieldDict[field_name]
        V = field.function_space()
        u = TrialFunction(V)
        v = TestFunction(V)
        f = Constant(0)

        u_n = Function(V)
        u_n.vector()[:] = field.vector()[:]
        F = u * v * dx + dt * alpha * dot(grad(u), grad(v)) * dx \
            - (u_n + dt * f) * v * dx
        a, L = lhs(F), rhs(F)

        # Time-stepping
        u = Function(V)
        t = 0
        bc = []
        vtkfile = File('heat_gaussian/solution.pvd')
        for n in range(n_timesteps):
            # Update current time
            t += dt

            # Compute solution
            solve(a == L, u, bc)

            vtkfile << (u, t)

            # Update previous solution
            u_n.assign(u)
        field.vector()[:] = u_n.vector()[:]
        self._fieldDict[field_name] = field

    def getVertexValues(self):
        """
        This function returns an instance of the class meanFlowVertexValues,
        which contains the vertex values of the fenics-Functions in the
        _fieldDict

        Function arguments:

        Function returns:
        instance of the class meanFlowVertexValues
        """
        return meanFlowVertexValues(self._meanfieldDict, self._oneFieldArray,self._FEMSpaces.exportMesh)

    def _getMeanFieldsToBeRead(self):
        listOfFieldsToBeRead = self._param.getMeanFlowFieldNames()
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadEnergy())
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadEoS())
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadReaction())
        
        # Delete duplicates
        listOfFieldsToBeRead = list(dict.fromkeys(listOfFieldsToBeRead))
        
        printDebug(True,"-- Mean flow fields to be read are: "+str(listOfFieldsToBeRead))
        return listOfFieldsToBeRead
    

class meanFlowVertexValues(fieldProperties):
    """
    This class provides the vertex values of the mean flow

    Parent classes:
    - fieldProperties

    Child classes:

    Private attributes:

    Protected attributes:

    - _fieldDict: The dictionary containing all the vertex values
    of the mean field

    Public attributes
    """

    def __init__(
            self, 
            fieldDict, 
            oneField,
            mesh,
            ):
        self._fieldDict = {}
        import numpy as np
        self._oneField = oneField.x.array[:]
        self._isMean = True
        for key in list(fieldDict.keys()):
            #tempMeanArray = fieldDict[key].compute_vertex_values()
            # The vector components (velocity u) need to be reshaped
            if fieldDict[key].function_space.num_sub_spaces > 1:
                tempSolutionArray = np.zeros((fieldDict[key].function_space.num_sub_spaces, mesh.coordinates().shape[0]),
                                             dtype=complex)
                for subSpace in range(fieldDict[key].function_space.num_sub_spaces):
                    indicesOfSubSpace = fieldDict[key].function_space.sub(subSpace).collapse()[1]
                    tempSolutionArray[subSpace, :] = fieldDict[key].x.array[indicesOfSubSpace]
                self._fieldDict[key] = tempSolutionArray
            else:
                self._fieldDict[key] = fieldDict[key].x.array[:]
