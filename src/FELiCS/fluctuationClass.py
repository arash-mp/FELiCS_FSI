# Standard libraries
from os import (
    listdir,
    remove,
)
from os.path import (
    exists
)

# Third party libraries
import numpy as np
from h5py import (
    File
)
from dolfinx.fem import (
    Function,
	FunctionSpace,
	Constant,
)
from ufl import (
    TrialFunctions,
)

# Local libraries and methods
from FELiCS.fieldProperties import fieldProperties
from FELiCS.dependentVariables.energyHandler import energyHandler
from FELiCS.dependentVariables.equationOfStateHandler import equationOfStateHandler
from FELiCS.dependentVariables.heatReleaseHandler import heatReleaseHandler
from FELiCS.dependentVariables.momentumHandler import momentumHandler
from FELiCS.dependentVariables.reactionHandler import reactionHandler
from FELiCS.functions import (
    printError,
    printWarning,
)
from FELiCS.export import export

from FELiCS.tensorUtils import (
    Tensor,
)

class fluctuationClass(
    fieldProperties,
    reactionHandler,
    equationOfStateHandler,
    energyHandler,
    momentumHandler,
    ):
    """
    This class fulfills two purposes. First, it is a wrapper for the test
    functions. These correspond to the actually transported variables in the
    equations and can be accessed from outside the class via public variables.
    In addition to these transported variables, this class also provides
    additional fluctuations in quantities, which are a function of the
    transported variables (dependent variables). These are also accessible via
    public variables. The second purpose is in the postprocessing. Since the
    linear algebra solution only consists of the transported variables, this
    class is used to determine the dependent variables from the transported
    variables, which can be used for further data processing/export.

    Parent classes:
    - fieldProperties
    - viscosityHandler
    - enthalpyHandler
    - heatReleaseHandler
    - laminarDiffusionHandler
    - temperatureHandler

    Child classes:

    Private attributes:
    - _fluc: The trial function in the mixed finite element space
    - _mean: The temporal mean flow (meanFlowClass)
    - _transportedQuantities: A List containing the string-names of all
        transported quantities
    - _zeroVelocityField: A fenics vector function containing zeros only

    Protected attributes:
    - _fieldDict: Dictionary containing both transported and dependent variables
    - _zeroField: A fenics function containing zeros only

    Public attributes:
    """

    def __init__(
            self,
            param,
            mean,
            FEMSpaces,
            coordinateSystem,
        ):
        """
        Constructor of the fluctuationClass. This function initializes the
        fields

        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object
        - FEMSpaces: FELiCS FEM spaces object
        - postProcessing: Flag if this class used in postprocessing or not

        Function returns:
        """
        self._param = param
        self._FEMSpaces = FEMSpaces
        self._coordinateSystem = coordinateSystem
        fieldProperties.__init__(
            self
        )
        self._isFluctuation = True
        self._isMean = False
        self._isSolution = False
        self._zeroField = Function(FEMSpaces.P2)
        self._zeroVelocityField \
            = Function(FEMSpaces.FunctionSpaceVectorVelocity)
        self._zeroField = Tensor(
                                       Function(self._FEMSpaces.P2),
                                       self._coordinateSystem,
                                       containsFluctuation = True,
                                       )
        self._fieldDict = {}
        self._mean = mean
        self._transportedQuantities = param.Case.getTransportedQuantityList()

        # _fluc is constructed. 
        self._fluc = TrialFunctions(FEMSpaces.VMixed)
        for field in self._transportedQuantities:
            indexOfFieldInList = self._transportedQuantities.index(field)
            self._fieldDict[field] = Tensor(
                                            self._fluc[indexOfFieldInList],
                                            self._coordinateSystem,
                                            containsFluctuation = True,
                                            )

        # Get all the variables, which need to be present
        neededVariables = []

        if not param.Case.SetOfEquations['Momentum']['Equation'] in ['None']:
            momentumHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearMomentum()

        if not param.Case.SetOfEquations['EquationOfState']['Equation'] in ['None']:
            equationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEoS()

        if not param.Case.SetOfEquations['Energy']['Equation'] in ['None']:
            energyHandler.__init__(
                                  self,
                                  )
            neededVariables += self._getNeededFieldsForLinearEnergy()
            

        # Delete duplicates
        neededVariables = list(dict.fromkeys(neededVariables))

        # While not all needed fluctuations are calculated, try calculating them
        n_try = 1
        while not set(neededVariables).issubset((self._fieldDict.keys())):
            if not param.Case.SetOfEquations['Momentum']['Equation'] in ['None']:
                self._relateConservativeToPrimitiveVariablesMomentum()
                self._initializeMolecularMomentumDiffusionFluctuation()
            if not param.Case.SetOfEquations['EquationOfState']['Equation'] in ['None']:
                self._initializeEoSFluctuations()
            if not param.Case.SetOfEquations['Energy']['Equation'] in ['None']:
                self._relateConservativeToPrimitiveVariablesEnergy()
                self._initializeMolecularHeatDiffusionFluctuation()
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                raise Exception('Attempt to calculate secondary variables not successful. Missing quantities: '\
                     + str(notInitializedFields))

        # print(self._param.Case.Mixture.getReactionMechanism()['type'])
        if not self._param.Case.Mixture.getReactionMechanism()['type'] == 'None':
            reactionHandler.__init__(
                self,
                )
            self._initializeReactions()
            
class fluctuationSolutions(
    fieldProperties,
    heatReleaseHandler,
    equationOfStateHandler,
    momentumHandler,
    energyHandler,
    reactionHandler,
    export,
):
    """
    This class contains the solutions to the linearized equations in form of the
    transported variables, which are stored in the vmmixed_vector.

    The class is can export and import the solutions in compact and readable
    format.

    TO DO: The export-method 'o', 'a' and 'a+' produce a high residuum at one
    vertex point, which lead to failure of the validation cases. The error
    leading to that residuum should be found

    Parent classes:
    - export
    - fieldProperties
    - viscosityHandler
    - enthalpyHandler
    - heatReleaseHandler
    - LaminarDiffusion
    - temperatureHandler

    Child classes:

    Private attributes:
    - _zeroVelocityField: A fenics vector function in the calculation space
        containing zeros only
    - _FEMSpaces: FEMSpaces Object. Contains the high-dimensional calculation
        space
    - _mean: meanflow object. Instance of meanflowclass
#	- _fieldDict: Contains the meanfield values in a dictionary
    - _transportedQuantities
    - _param: FELiCS parameter object
    - _meshfilename: filename of the h5-file including the exportMesh
    - _linearFunctionSpaces: FEMSpaces Object, contains the FEMSpaces of first
        order of the exportMesh
    - _exportZeroScalarField: fenics function containing a zeroField on the
        exportMesh
    - _exportZeroVectorField: fenics function containing a zero vector Field on
        the exportMesh
    - _exportMesh: Object of type FELiCSMesh, contains the export Mesh.
    - _uValuesList: list of velocity component names
    - _solution: a list containing all attributes of the solution
    - _meanfieldDict: dictionary containing Field-Functions, which are
        interpolated to the export Space.

    Protected attributes:
    - _fieldDict: Dictionary containing both transported and dependent variables
    - _zeroField: A fenics function containing zeros only

    """

    def __init__(
            self,
            param,
            mean,
            FEMSpaces,
            omega,
            vmixedVector,
            isResponseOrDirect,
            gainNumber=-1,
            gainValue=-1
    ):
        """
        Function arguments:
        - param: FELiCS parameter object
        - mean: FELiCS mean flow object
        - FEMSpaces: FELiCS FEM spaces object
        - omega: complex Eigenvalue of the solution
        - vmixedVector: numpy-array containing the complex values of the
            transported quantities at every DOF-coordinate.
        - isResponseOrDirect: Boolean, is true if the solution is of Kind
            Response or Direct and false if it is of type Forcing or Adjoint.
        - gainNumber: Integer, which contains the number of the gain. For
            solutions, which have no gain it is set -1

        Function returns:
        """

        # dolfinx specific: There is no compute_vertex_values anymore.
        self._isSolution = True
        self._zeroField = Function(FEMSpaces.P2)
        self._zeroField.x.array[:] = 0.0
        self._zeroField = self._zeroField.x.array[:]
        #self._zeroVelocityField = Function(FEMSpaces.FunctionSpaceVectorVelocityP1)
        #self._zeroVelocityField.x.array[:] = 0.0

        # if np.imag(gainValue) > 1e-10 * np.real(gainValue):
        #     printWarning('The gain is a complex number, while it should be \
        #     real. I will ignore this and take the real part!')
        self._gainValue = np.real(gainValue)
        self._omega = omega
        self._isResponseOrDirect = isResponseOrDirect
        self._vmixedVector = vmixedVector
        self._gainNumber = gainNumber
        self._FEMSpaces = FEMSpaces
        self._mean = mean
        # self._fieldDict = mean.fieldDict
        self._transportedQuantities = param.Case.getTransportedQuantityList()
        self._param = param
        export.__init__(
            self,
            param,
            FEMSpaces)

    def _flucExportWrapper(
            self,
            group):
        """
        This function calls the methods for the fluc export. Therefore, it
        iterates over the sub-solutions and creates subgroups if needed.

        Function arguments:
        - group: Points on a group inside the h5-file, in which the solution
            should be exported.

        Function returns:

        """

        self._fieldDict = self._mapCalcToExport(self._vmixedVector)

        # Get all the variables, which need to be present
        neededVariables = []
        if not self._param.Case.SetOfEquations['Momentum']['Equation'] in ['None']:
            momentumHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearMomentum()

        if not self._param.Case.SetOfEquations['EquationOfState']['Equation'] in ['None']:
            equationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEoS()

        if not self._param.Case.SetOfEquations['Energy']['Equation'] in ['None']:
            energyHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEnergy()
    
        if not self._param.Case.SetOfEquations['Energy']['Equation'] in ['None']:
            energyHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEnergy()
    
        # Delete duplicates
        neededVariables = list(dict.fromkeys(neededVariables))

        # While not all needed fluctuations are calculated, try calculating them
        n_try = 0
        while not set(neededVariables).issubset((self._fieldDict.keys())):
            meanVertexValues = self._mean.getVertexValues()
            if not self._param.Case.SetOfEquations['Momentum']['Equation'] in ['None']:
                self._relateConservativeToPrimitiveVariablesMomentum(
                                                    meanVertexValues,
                                                    )
                self._initializeMolecularMomentumDiffusionFluctuation(meanVertexValues)

            if not self._param.Case.SetOfEquations['EquationOfState']['Equation'] in ['None']:
                self._initializeEoSFluctuations(
                                            meanVertexValues,
                                                )

            if not self._param.Case.SetOfEquations['Energy']['Equation'] in ['None']:
                self._relateConservativeToPrimitiveVariablesEnergy(
                                            meanVertexValues,
                                                )
                self._initializeMolecularHeatDiffusionFluctuation(
                                            meanVertexValues,
                                            )
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                raise Exception('Attempt to calculate secondary variables not successful. Missing quantities: ' +notInitializedFields)

        if not self._param.Case.Mixture.getReactionMechanism()['type'] == 'None':
            reactionHandler.__init__(
                self,
                )
            self._initializeReactions(
                meanVertexValues,
                )
        self._writeDictToH5(
            self._fieldDict,
            group,
            True
        )

    def _importSolVector(
            self,
            filename
    ):
        """
        This function imports the VMixed-solution vector if it is present in raw
        form inside a h5-file with the name filename.

        Function arguments:
        - filename: String filename of the h5-file. To mark it as including the
            VMixed-Vector, it ends with "_sol.h5"

        Function returns:
        - returns complex valued numpy array, which represents the
            VMixed-solution

        """
        hf = File(filename, 'r')
        frequency = hf['fluctuation/0/'].attrs.get('frequency')
        fieldMagnitude = np.array(hf[f'fluctuation/0/{frequency}/magnitude'][:])
        fieldAngle = np.array(hf[f'fluctuation/0/{frequency}/angle'][:])

        return fieldMagnitude * np.exp(1j * fieldAngle)

    def exportSolution(
            self,
            filename,
            flag
    ):
        """
        This method exports the Solution to the hdf5-format.
        Furthermore, a XMF-file is constructed with which the solution can be
        opened in paraview.

        Function arguments:
        - filename: filename of the h5-file, where the fields should be exported
        - flag: string flag, specifying the kind of export. Possible values are
            'c', 'o', 'a' and 'a+'

        Function returns:

        """
        filenameWithoutExtension = filename.split('.h5')[0]
        filenameWithoutFolder = filename
        filename = f'{self._param.Export.ExportFolder}/' + filename

        # appendFlag = False

        if flag == 'o':

            if exists(filename):
                remove(filename)
            hf = File(filename, 'w')
            pointGroup = self._createH5GroupStructure(hf)

            self._flucExportWrapper(pointGroup)

            self.writeXMFFile(f'{self._param.Export.ExportFolder}/'
                              + f'{self._param.Case.AnalysisMode}_mesh.h5',
                              self._mean.meanflowFilename, hf)

            # export the param-object to the h5-file as string:
            self._param.export(f'{hf.filename}')

            hf.close()

        elif flag == 'a':

            # open the file in append mode:
            hf = File(filename, 'a')
            pointGroup = self._createH5GroupStructure(hf)

            # appendFlag = True

            self._flucExportWrapper(pointGroup)

            self.writeXMFFile(f'{self._param.Export.ExportFolder}/'
                              + f'{self._param.Case.AnalysisMode}_mesh.h5',
                              self._mean.meanflowFilename, hf)

            # export the param-object to the h5-file as string:
            self._param.export(f'{hf.filename}')

            hf.close()
        elif flag == 'a+':
            i = 0
            fileNameWithNumb = filenameWithoutExtension + '_0' + '.h5'
            fileNameListing = listdir(self._param.Export.ExportFolder)

            while fileNameWithNumb in fileNameListing:
                i += 1
                fileNameWithNumb = filenameWithoutExtension + '_' + f'{i}' \
                                   + '.h5'

            nextFileindex = i
            hf = File(f'{self._param.Export.ExportFolder}/'
                      + filenameWithoutExtension + '_' + f'{i}' + '.h5', 'w')

            pointGroup = self._createH5GroupStructure(hf)

            self._flucExportWrapper(pointGroup)

            self.writeXMFFile(f'{self._param.Export.ExportFolder}/'
                              + f'{self._param.Case.AnalysisMode}_mesh.h5',
                              self._mean.meanflowFilename, hf)

            # export the param-object to the h5-file as string:
            self._param.export(f'{hf.filename}')

            hf.close()
        elif flag == 'c':
            # writes out the VMixed-Vector to File:

            # add suffix 'sol' to filename:
            filenameForCExport = filenameWithoutExtension + '_sol.h5'

            # if the file already exists, overwrite it:
            hf = File(f'{self._param.Export.ExportFolder}/{filenameForCExport}',
                      'w')

            pointGroup = self._createH5GroupStructure(hf)
            frequency = pointGroup.parent.attrs.get('frequency')
            frequencyGroup = pointGroup.parent[frequency]

            frequencyGroup.create_dataset('magnitude',
                                          data=np.abs(self._vmixedVector))
            frequencyGroup.create_dataset('angle',
                                          data=np.angle(self._vmixedVector))
            frequencyGroup.file.close()

            # export the param-object to the h5-file as string:
            self._param.export(f'{self._param.Export.ExportFolder}\
            /{filenameForCExport}')

        else:
            print('The given flag is not known!')

    def importSolution(
            self,
            filename
    ):
        """
        This Function imports a solution into the VMixed space.

        Function arguments:
        - filename: String containing the filename of the h5-File.

        Function Returns:

        """
        # if the filename has the suffix "_sol", its the raw VMixed-Vector and
        # the import is easy:
        if '_sol.h5' in filename:
            importedSolVector = self._importSolVector(
                f'{self._param.Export.ExportFolder}/{filename}')
            self._param.importFromFile(
                f'{self._param.Export.ExportFolder}/{filename}')
            assert np.allclose(importedSolVector, self._vmixedVector)
            self._vmixedVector = importedSolVector

    @property
    def solutVector(self):
        return self._vmixedVector

    @property
    def solutionKind(self):
        if self._param.Case.AnalysisMode == 'Resolvent':
            if self._isResponseOrDirect:
                return 'Response'
            else:
                return 'Forcing'
        elif self._param.Case.AnalysisMode == 'Input-Output':
            if self._isResponseOrDirect:
                return 'Response'
            else:
                return 'Forcing'
        elif self._param.Case.AnalysisMode == 'Modal':
            if self._isResponseOrDirect:
                return 'Direct'
            else:
                return 'Adjoint'
        else:
            return ""

    @property
    def omega(self):
        return self._omega

    @property
    def gainNumber(self):
        return self._gainNumber

    @property
    def gainValue(self):
        return self._gainValue
