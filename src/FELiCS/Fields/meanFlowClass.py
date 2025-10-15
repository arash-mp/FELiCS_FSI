# Third party libraries
import h5py
import numpy as np

# Local libraries and methods
from FELiCS.IO.export                                           import export
from FELiCS.Fields.fieldProperties                              import fieldProperties
from FELiCS.Equation.dependentVariables.energyHandler           import energyHandler
from FELiCS.Equation.dependentVariables.equationOfStateHandler  import equationOfStateHandler
from FELiCS.Equation.dependentVariables.reactionHandler         import reactionHandler
from FELiCS.Misc.logging                                        import Logger
from FELiCS.IO.Reader                                           import Reader
from FELiCS.IO.Mapping                                          import Mapping
from FELiCS.Fields.Field                                        import Field

# Get the logger
logger = Logger.get_logger("felics")

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
    def __init__(self, param, FEMSpaces, mesh):
        
        # Initialization
        self._isMean            = True
        self._isFluctuation     = False
        self._param             = param
        self._FEMSpaces         = FEMSpaces
        self._mesh              = mesh
        self._coordinateSystem  = mesh.coordinateSystem
        self._meanflowFilename  = None
        self._mixture           = param.Mixture
        
        # Define some useful tensors
        self._zeroVectorField   = Field(self._FEMSpaces.FunctionSpaceVectorVelocity, self._mesh, name="zero")
        self._zeroField         = Field(self._FEMSpaces.P2, self._mesh, name="zero")
        self._oneField          = Field(self._FEMSpaces.P2, self._mesh, name="one")
        self._oneField.setConstantValue(1.0)

        # Allocate
        self._customMeanFlowQuantities = []
        

    def importDataFromFileAndExportToH5(self):
        logger.info(f"Reading input flow from: '{self._param.FlowInput.MeanFlowFilePath}'")

        # Initialization
        self._fieldDict             = {}
        self._notInFileList         = []
        self._ScalarFunctionSpace   = self._FEMSpaces.P2
        self._VectorFunctionSpace   = self._FEMSpaces.FunctionSpaceVectorVelocity
        
        # Get variable list to be read
        nameListMean                = self._getMeanFieldsToBeRead()
        
        # Get the group name depending on the file type
        # TODO: remove groups into FELiCS files
        if self._param.Case.needInterpolation:
            groupName               = "MeanFlow"
        else:
            groupName               = "meanflow"

        # Initialize the reader
        reader          = Reader(
            needInterpolation   = self._param.Case.needInterpolation,
            originalMeshFile    = None,
            isComplex           = False, # NOTE: Should be an attribute of the field, not the reader?
            cacheData           = True,
        )
        # TEMPORARY: Bind the params and FEM spaces to the reader
        reader.bind_env(self._param, self._FEMSpaces)

        # Create empty fields for each variable 
        # TODO: Create new type of Field collection for mean flow? (Similar to mode collection.)
        for name in nameListMean:
            # TODO: Include the tensor order in the field names to avoid hardcoding.
            if name == 'u' or name == 'u_forcing_r' or name == 'u_forcing_i':
                field = Field(
                    self._FEMSpaces.FunctionSpaceVectorVelocity,
                    self._mesh,
                    name
                )
            else:
                field = Field(
                    self._FEMSpaces.P2,
                    self._mesh,
                    name
                )
                
            # Load the data from file into the field
            self._fieldDict[name], notInFile = field.importData(
                reader,
                self._param.FlowInput.MeanFlowFilePath,
                groupName
            )
            
            # Variables not found in the file are stored in a list
            if notInFile:
                self._notInFileList.append(name)
        
        # set 'ut' field to zero
        if 'ut' in list(self._fieldDict.keys()):
            self._fieldDict['ut'].setConstant(0.)
                
        # Define the viscosity and alfa fields
        # NOTE: This should move to a handler
        self.initLamDiff()
        self.initThermodynamicQuantities()

        # export mean flow to "meanflow.h5" file
        self.mapToExportMeshAndExport(self._FEMSpaces, "MeanFlow.h5")


    def initLamDiff(self):
        if self._param.Case.MolViscModel == 'Constant':
            self._fieldDict['nulam']            = Field(self._FEMSpaces.P2, self._mesh, name = "nulam")
            self._fieldDict['nulam'].setCoefficientArray(self._param.Case.MolVisc)

        for specie in self._param.Mixture.getSpeciesList('transported'):
            Sc      = self._param.Mixture.species[specie]['Sc']
            nuTot   = Field(self._ScalarFunctionSpace, self._mesh, name = "nuTot")

            if 'nulam' in list(self._fieldDict.keys()):
                nuTot += self._fieldDict['nulam']
            if 'nuturb' in list(self._fieldDict.keys()):
                nuTot += self._fieldDict['nuturb']
            if 'nuSGS' in list(self._fieldDict.keys()):
                nuTot += self._fieldDict['nuSGS']
                
            self._fieldDict['D_' + specie]  = Field(self._FEMSpaces.P2, self._mesh, name= 'D_'+specie)
            self._fieldDict['D_' + specie]  = nuTot / Sc

    def initThermodynamicQuantities(self):
        self._fieldDict['Pr']            = Field(self._FEMSpaces.P2, self._mesh, name="Pr")
        self._fieldDict['Pr'].setConstantValue(self._param.Case.PrandtlNumber)

    def mapToExportMeshAndExport(self, FEMSpaces, filename):
        logger.info("Mapping mean flow to export mesh and exporting.")
        self._meanflowFilename  = filename
        filehandler             = h5py.File(f'{self._param.Export.ExportFolder}/{filename}', 'w')
        group                   = filehandler.create_group('meanflow')
        export.__init__(self, self._param, self._FEMSpaces)
        #self._meanfieldDict, dictImag = self._mapCalcToExport(self._fieldDict)
        self._meanfieldDict     = self._mapCalcToExport(self._fieldDict)
        #exportDict = self._calculateVertexValuesFromDict(self._meanfieldDict,
        #                                                 dictImag)
        #self._writeDictToH5(exportDict, group)
        self._writeDictToH5(
            self.getVertexValues()._fieldDict,
            group,
        )
        filehandler.close()

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

        self._hSpec = {}
        self._hSpec['CH4'] = project(janaf.janaf_Hs_expr(CH4, self.T),
                                      self._FEMSpaces.P2)
        self._hSpec['O2'] = project(janaf.janaf_Hs_expr(O2, self.T),
                                     self._FEMSpaces.P2)
        self._hSpec['CO'] = project(janaf.janaf_Hs_expr(CO, self.T),
                                     self._FEMSpaces.P2)
        self._hSpec['CO2'] = project(janaf.janaf_Hs_expr(CO2, self.T),
                                      self._FEMSpaces.P2)
        self._hSpec['H2O'] = project(janaf.janaf_Hs_expr(H2O, self.T),
                                      self._FEMSpaces.P2)

    def getVertexValues(self):
        """
        This function returns an instance of the class meanFlowVertexValues,
        which contains the vertex values of the fenics-Functions in the
        _fieldDict

        Function arguments:

        Function returns:
        instance of the class meanFlowVertexValues
        """
        return meanFlowVertexValues(self._meanfieldDict, self._FEMSpaces.exportMesh)

    def _getMeanFieldsToBeRead(self):
        """
        Generate a list of mean flow field names to be read from the data source.
        This method combines field names from various sources:
        - Basic mean flow fields from parameters
        - Additional fields needed for energy calculations
        - Additional fields needed for equation of state calculations
        - Additional fields needed for reaction calculations
        - Custom mean flow quantities defined by the user
        It also removes any duplicate field names to ensure each field is only read once.
        TODO: Include the tensor order in the field names
        Returns:
            list: A deduplicated list of field names to be read
        """
        listOfFieldsToBeRead = self._param.getMeanFlowFieldNames()
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadEnergy())
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadEoS())
        listOfFieldsToBeRead.extend(self._additionalFieldsToBeReadReaction())
       
        listOfFieldsToBeRead.extend(self._customMeanFlowQuantities[:])

        # Delete duplicates
        listOfFieldsToBeRead = list(dict.fromkeys(listOfFieldsToBeRead))
        
        logger.debug("Mean flow fields to read: "+str(listOfFieldsToBeRead))
        return listOfFieldsToBeRead

    def addCustomMeanFlowQuantity(self,key):
        self._customMeanFlowQuantities.append(key)

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

    def __init__(self,fieldDict,mesh):
        
        self._fieldDict     = {}
        self._isMean        = True
        
        # Define some useful tensors
        self._zeroField     = np.zeros(mesh.coordinates().shape[0])
        self._oneField      = np.ones(mesh.coordinates().shape[0])
        
        for key in list(fieldDict.keys()):
            #tempMeanArray = fieldDict[key].compute_vertex_values()
            # The vector components (velocity u) need to be reshaped
            if fieldDict[key].space.num_sub_spaces > 1:
                tempSolutionArray = np.zeros((fieldDict[key].space.num_sub_spaces, mesh.coordinates().shape[0]),
                                             dtype=complex)
                for subSpace in range(fieldDict[key].space.num_sub_spaces):
                    indicesOfSubSpace = fieldDict[key].space.sub(subSpace).collapse()[1]
                    tempSolutionArray[subSpace, :] = fieldDict[key].function.x.array[indicesOfSubSpace]
                self._fieldDict[key] = tempSolutionArray
            else:
                self._fieldDict[key] = fieldDict[key].getCoefficientArray()
