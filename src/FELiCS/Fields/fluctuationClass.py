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
)
from ufl import (
    TrialFunctions,
)

# Local libraries and methods
from    FELiCS.Fields.fieldProperties                               import fieldProperties
from    FELiCS.Fields.Field                                         import Field
from    FELiCS.Equation.dependentVariables.energyHandler            import energyHandler
from    FELiCS.Equation.dependentVariables.equationOfStateHandler   import equationOfStateHandler
from    FELiCS.Equation.dependentVariables.heatReleaseHandler       import heatReleaseHandler
from    FELiCS.Equation.dependentVariables.momentumHandler          import momentumHandler
from    FELiCS.Equation.dependentVariables.reactionHandler          import reactionHandler
from    FELiCS.IO.export                                            import export
from    FELiCS.Misc.tensorUtils                                     import Tensor
from 	FELiCS.Misc.logging                                         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class fluctuationClass(
    fieldProperties,
    reactionHandler,
    equationOfStateHandler,
    energyHandler,
    momentumHandler,
    ):
    """
    Container for fluctuating variables and derived fields for postprocessing.

    This class wraps the test functions representing transported variables
    within the linearized governing equations. In addition to these, it also
    derives secondary fields (dependent variables) from the transported
    variables, which are needed for postprocessing and exporting results.

    It plays a central role in calculating and organizing both transported
    and derived quantities used in numerical simulations.

    **Initialize the fluctuationClass object**

    Parameters
    ----------
    param : FELiCSParameter
        FELiCS parameter object containing simulation configuration.
    mean : meanFlowClass
        Mean flow object.
    FEMSpaces : FEMSpaceHandler
        Object encapsulating FEM spaces used in the simulation.
    coordinateSystem : object
        Representation of the simulation's coordinate system.

    Attributes
    ----------
    _fluc : ufl.argument.TrialFunction
        Trial function for the mixed finite element space.
    _mean : meanFlowClass
        Temporal mean flow object.
    _transportedQuantities : list of str
        Names of the transported quantities.
    _zeroVectorField : dolfinx.Function
        Zero-valued vector function in the velocity space.
    _fieldDict : dict
        Dictionary of calculated fields (both transported and dependent).
    _zeroField : dolfinx.Function
        Scalar zero field in the scalar space.
    """

    def __init__(
            self,
            param,
            mean,
            FEMSpaces,
            coordinateSystem,
        ):
        """
        Initializes fluctuation fields and derives dependent variables.

        Parameters
        ----------
        param : FELiCSParameter
            FELiCS parameter object.
        mean : meanFlowClass
            Mean flow object.
        FEMSpaces : FEMSpaceHandler
            FEM spaces to be used in calculation.
        coordinateSystem : object
            Coordinate system used in the domain.
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
        # TODO Sophie: set correct name
        self._zeroVectorField \
            = Field(FEMSpaces.FunctionSpaceVectorVelocity, self._param.getMesh(), name=[]).getTensor()
        self._zeroField = Field(FEMSpaces.P2, self._param.getMesh(), name="zero").getTensor()
        self._fieldDict = {}
        self._mean = mean
        self._transportedQuantities = param.getTransportedQuantityList()

        # _fluc is constructed. 
        self._fluc = TrialFunctions(FEMSpaces.VMixed)
        for field in self._transportedQuantities:
            indexOfFieldInList = self._transportedQuantities.index(field)
            self._fieldDict[field] = Tensor(
                                            self._fluc[indexOfFieldInList],
                                            self._coordinateSystem,
                                            mayHaveSpectralDimension = True,
                                            )

        # Get all the variables, which need to be present
        neededVariables = []

        if param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
            momentumHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearMomentum()

        if param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
            equationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEoS()

        if param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
            energyHandler.__init__(
                                  self,
                                  )
            neededVariables += self._getNeededFieldsForLinearEnergy()
            

        # Delete duplicates
        neededVariables = list(dict.fromkeys(neededVariables))

        # While not all needed fluctuations are calculated, try calculating them
        n_try = 1
        while not set(neededVariables).issubset((self._fieldDict.keys())):
            if param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
                self._relateConservativeToPrimitiveVariablesMomentum()
                self._initializeMolecularMomentumDiffusionFluctuation()
            if param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
                self._initializeEoSFluctuations()
            if param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
                self._relateConservativeToPrimitiveVariablesEnergy()
                self._initializeMolecularHeatDiffusionFluctuation()
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                logger.error('Attempt to calculate secondary variables not successful. Missing quantities: '\
                     + str(notInitializedFields) + ". Maybe the mixture file is still in the old format (ending with a '.mix' instead of '.json')?")

        # print(self._param.Mixture.getReactionMechanism()['type'])
        if not self._param.Mixture.getReactionMechanism()['type'] == 'None':
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
    Stores and manages linearized fluctuation solutions and their export.

    This class contains the computed fluctuation solutions in the mixed
    function space. It supports exporting to and importing from HDF5/XDMF
    formats for postprocessing and visualization in tools like ParaView.

    **Initialize the fluctuationSolutions object**

    Parameters
    ----------
    param : FELiCSParameter
        FELiCS parameter object containing simulation configuration.
    mean : meanFlowClass
        Mean flow object.
    FEMSpaces : FEMSpaceHandler
        Finite element space handler.
    omega : complex
        Complex eigenvalue representing the frequency of the solution.
    vmixedVector : np.ndarray
        Complex-valued vector representing the mixed solution.
    isResponseOrDirect : bool
        Flag indicating whether the solution is a Response/Direct (True)
        or Forcing/Adjoint (False).
    gainNumber : int, optional
        Index of the gain value, default is -1.
    gainValue : float, optional
        Value of the gain, default is -1.

    Attributes
    ----------
    _zeroVectorField : dolfinx.Function
        Vector-valued zero field used for initialization (not stored).
    _FEMSpaces : FEMSpaceHandler
        High-order FEM spaces for simulation.
    _mean : meanFlowClass
        Mean flow object.
    _transportedQuantities : list of str
        Names of the transported quantities.
    _param : FELiCSParameter
        Simulation configuration object.
    _meshfilename : str
        Filename of the export mesh (if applicable).
    _linearFunctionSpaces : FEMSpaceHandler
        Low-order FEM spaces for export mesh (if applicable).
    _exportZeroScalarField : dolfinx.Function
        Scalar zero field on export mesh.
    _exportZeroVectorField : dolfinx.Function
        Vector zero field on export mesh.
    _exportMesh : FELiCSMesh
        Export mesh.
    _uValuesList : list of str
        List of velocity component names.
    _solution : list
        Solution vector components.
    _meanfieldDict : dict
        Dictionary containing mean field interpolated onto export space.
    _fieldDict : dict
        Dictionary of calculated fields (both transported and dependent).
    _zeroField : np.ndarray
        Zero-valued array for initialization.
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
        Initialize the fluctuationSolutions object.

        Parameters
        ----------
        param : FELiCSParameter
            FELiCS parameter object.
        mean : meanFlowClass
            Mean flow object.
        FEMSpaces : FEMSpaceHandler
            Finite element space handler.
        omega : complex
            Complex eigenvalue of the solution.
        vmixedVector : np.ndarray
            Complex-valued vector representing the mixed solution.
        isResponseOrDirect : bool
            Flag indicating whether the solution is of type Response/Direct or Forcing/Adjoint.
        gainNumber : int, optional
            Index for gain tracking, default is -1.
        gainValue : float, optional
            Gain value associated with the solution, default is -1.
        """

        # dolfinx specific: There is no compute_vertex_values anymore.
        self._isSolution = True
        self._zeroField = Function(FEMSpaces.P2)
        self._zeroField.x.array[:] = 0.0
        self._zeroField = self._zeroField.x.array[:]
        #self._zeroVectorField = Function(FEMSpaces.FunctionSpaceVectorVelocityP1)
        #self._zeroVectorField.x.array[:] = 0.0

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
        self._transportedQuantities = param.getTransportedQuantityList()
        self._param = param
        export.__init__(
            self,
            param,
            FEMSpaces)

    def _flucExportWrapper(self, group):
        """
        Internal method to compute and export all fluctuation fields.

        This method maps the solution vector to field data, derives
        dependent variables if needed, and writes everything to the
        provided HDF5 group.

        Parameters
        ----------
        group : h5py.Group
            HDF5 group in which fields are to be written.
        """

        self._fieldDict = self._mapCalcToExport(self._vmixedVector)

        # Get all the variables, which need to be present
        neededVariables = []
        if self._param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
            momentumHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearMomentum()

        if self._param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
            equationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEoS()

        if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
            energyHandler.__init__(
                self,
                )
            neededVariables += self._getNeededFieldsForLinearEnergy()
    
        if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
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
            if self._param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
                self._relateConservativeToPrimitiveVariablesMomentum(
                                                    meanVertexValues,
                                                    )
                self._initializeMolecularMomentumDiffusionFluctuation(meanVertexValues)

            if self._param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
                self._initializeEoSFluctuations(
                                            meanVertexValues,
                                                )

            if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
                self._relateConservativeToPrimitiveVariablesEnergy(
                                            meanVertexValues,
                                                )
                self._initializeMolecularHeatDiffusionFluctuation(
                                            meanVertexValues,
                                            )
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                logger.error('Attempt to calculate secondary variables not successful. Missing quantities: ' +notInitializedFields +". Maybe the mixture file is still in the old format (ending with a '.mix' instead of '.json')?")

        if not self._param.Mixture.getReactionMechanism()['type'] == 'None':
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

    def _importSolVector(self, filename):
        """
        Import raw VMixed solution vector from HDF5 file.

        Parameters
        ----------
        filename : str
            Path to the HDF5 file ending in "_sol.h5".

        Returns
        -------
        np.ndarray
            Complex-valued numpy array representing the solution.
        """
 
        hf = File(filename, 'r')
        frequency = hf['fluctuation/0/'].attrs.get('frequency')
        fieldMagnitude = np.array(hf[f'fluctuation/0/{frequency}/magnitude'][:])
        fieldAngle = np.array(hf[f'fluctuation/0/{frequency}/angle'][:])

        return fieldMagnitude * np.exp(1j * fieldAngle)

    def exportSolution(self, filename, flag):
        """
        Export fluctuation solution to HDF5/XMF format.

        Parameters
        ----------
        filename : str
            Base filename for the exported HDF5 file.
        flag : str
            Export mode. Options are:
            - 'o': overwrite
            - 'a': append
            - 'a+': write with new index suffix
            - 'c': compact, raw vector only

        Raises
        ------
        ValueError
            If the provided flag is unknown.
        """

        filenameWithoutExtension = filename.split('.h5')[0]
        # filenameWithoutFolder = filename
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

            # nextFileindex = i
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

    def importSolution(self, filename):
        """
        Import a previously exported VMixed fluctuation solution.

        Parameters
        ----------
        filename : str
            Filename of the HDF5 file to import from.
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
        """
        Complex-valued vector of the fluctuation solution.

        Returns
        -------
        np.ndarray
            Solution vector in mixed function space.
        """
        return self._vmixedVector


    @property
    def solutionKind(self):
        """
        Type of the solution based on analysis mode and configuration.

        Returns
        -------
        str
            One of 'Response', 'Forcing', 'Direct', 'Adjoint', or empty string.
        """

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
        """
        Eigenvalue (frequency) of the fluctuation solution.

        Returns
        -------
        complex
            Complex frequency associated with the solution.
        """
        return self._omega


    @property
    def gainNumber(self):
        """
        Gain number identifier.

        Returns
        -------
        int
            Gain index or -1 if unused.
        """
        return self._gainNumber


    @property
    def gainValue(self):
        """
        Gain magnitude.

        Returns
        -------
        float
            Real-valued gain.
        """
        return self._gainValue

