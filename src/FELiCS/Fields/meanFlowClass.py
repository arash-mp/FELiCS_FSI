#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# Third party libraries
import os
import h5py
import numpy as np

# Local libraries and methods
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
):
    """
    Container and handler for mean-flow fields on a given mesh.

    This class manages reading, storing, and post-processing of mean-flow
    quantities from external files and provides auxiliary fields and
    thermodynamic/transport properties required by FELiCS solvers.

    **Initialize the meanFlowClass object**

    Parameters
    ----------
    param : object
        Global parameter object providing case, mixture and I/O configuration,
        e.g. `param.FlowInput`, `param.Case`, and `param.Mixture`.
    FEMSpaces : object
        Collection of finite element function spaces, expected to provide
        attributes such as ``P2`` and ``FunctionSpaceVectorVelocity``.
    mesh : object
        Mesh object used to define the function spaces, expected to provide
        attributes such as ``coordinateSystem`` and ``exportMesh``.

    Attributes
    ----------
    _isMean : bool
        Flag indicating that this instance represents a mean flow (always True).
    _isFluctuation : bool
        Flag indicating that this instance represents a fluctuation field
        (always False for mean flow).
    _param : object
        Reference to the global parameter object used for configuration.
    _FEMSpaces : object
        Reference to the FEM spaces used for field definitions.
    _mesh : object
        Reference to the mesh used for the fields.
    _coordinateSystem : str or object
        Coordinate system information taken from the mesh.
    _meanflowFilename : str or None
        Optional path to the processed mean-flow file (set when exporting).
    _mixture : object
        Mixture model taken from ``param.Mixture``.
    _zeroVectorField : Field
        Vector field initialized to zero on the velocity space.
    _zeroField : Field
        Scalar field initialized to zero on the scalar space.
    _oneField : Field
        Scalar field initialized to one on the scalar space.
    _customMeanFlowQuantities : list of str
        List of additional mean-flow quantity names requested by the user.
    _fieldDict : dict of {str: Field}
        Dictionary of all mean-flow fields after import.
    _notInFileList : list of str
        Names of requested fields that were not found in the input file.
    _ScalarFunctionSpace : object
        Scalar function space used for viscosity and diffusion fields.
    _VectorFunctionSpace : object
        Vector function space used for velocity-related fields.
    xmfHeader : object
        Metadata returned by the writer when exporting mean fields.

    Notes
    -----
    The class delegates the specification of additional required fields to the
    handler base classes (`energyHandler`, `equationOfStateHandler`,
    `reactionHandler`). These base classes provide lists of extra field names
    needed for their respective models.
    """
    def __init__(self, param, FEMSpaces, mesh):
        """
        Initialize the meanFlowClass instance.

        Parameters
        ----------
        param : object
            Global parameter object providing case, mixture and I/O settings.
        FEMSpaces : object
            Object holding the finite element spaces used to construct fields.
        mesh : object
            Mesh object underlying all finite element spaces.
        """
        
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
        

    def importDataFromFileAndExportToH5(self, writer):
        """
        Import mean-flow fields from the configured input file and export them.

        This method reads the mean-flow quantities from the file specified in
        ``param.FlowInput.MeanFlowFilePath``, constructs `Field` instances for
        all required variables, computes additional transport and thermodynamic
        quantities, and exports the assembled set of fields to an HDF5 file.

        Parameters
        ----------
        writer : object
            Writer object providing an ``exportListOfFieldsToH5`` method used
            to export the list of mean-flow fields.

        Notes
        -----
        The set of fields to be read is assembled from:
        
        - the mean-flow field names defined in ``param``,
        - additional fields required by the energy, equation-of-state, and reaction handlers and 
        - any custom mean-flow quantities added by the user.

        If the density field ``rho`` is requested but not present in the file,
        its values are initialized to a constant of 1.0. If a turbulent
        velocity field ``ut`` exists, it is initialized to zero.

        After import, transport coefficients and thermodynamic quantities are
        initialized via :meth:`initLamDiff` and
        :meth:`initThermodynamicQuantities`.
        """
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
        reader                      = Reader(
            sourceDir               = os.path.dirname(self._param.FlowInput.MeanFlowFilePath),
            needInterpolation       = self._param.Case.needInterpolation,
            felicsMeshFilePath      = None if self._param.Case.needInterpolation else f'{self._param.Export.ExportFolder}/mesh.h5',
            cacheData               = True,
        )

        # Create empty fields for each variable 
        # TODO: Create new type of Field collection for mean flow? (Similar to mode collection.)
        for name in nameListMean:
            # TODO: Include the tensor order in the field names to avoid hardcoding.
            if name == 'u' or name == 'u_forcing_r' or name == 'u_forcing_i' or name == 'u_target':
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

            # correct value of rho, if not given, to 1
            if name == "rho" and np.linalg.norm(field.getCoefficientArray()) == 0:
                self._fieldDict[name].setConstantValue(1.)
            
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

        # export a list of all mean fields and export them in a "MeanFlow" file
        exportFields = []
        for key, value in self._fieldDict.items():
            exportFields.append(value)

        self.xmfHeader = writer.exportListOfFieldsToH5(exportFields, "MeanFlow")



    def initLamDiff(self):
        """
        Initialize molecular viscosity and species diffusion coefficients.

        Based on the selected molecular viscosity model and the mixture data,
        this method constructs the laminar viscosity field ``nulam`` (if a
        constant model is used) and computes species diffusion coefficients
        ``D_<specie>`` using the corresponding Schmidt numbers.

        The effective kinematic viscosity for each species is assembled from
        the available contributions:
        ``nulam``, ``nuturb`` and ``nuSGS``, if present in ``_fieldDict``.

        Notes
        -----
        This routine populates the following entries in ``_fieldDict``:

        - ``'nulam'`` : Field
            Constant laminar viscosity (for constant viscosity models).
        - ``'D_<specie>'`` : Field
            Species diffusion coefficients for transported species.
        """
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
        """
        Initialize thermodynamic mean-flow quantities.

        Currently this method creates and initializes the Prandtl number field
        ``'Pr'`` in ``_fieldDict`` using the value from
        ``param.Case.PrandtlNumber``.

        Notes
        -----
        The ``'Pr'`` field is defined on the scalar P2 function space.
        """
        self._fieldDict['Pr']            = Field(self._FEMSpaces.P2, self._mesh, name="Pr")
        self._fieldDict['Pr'].setConstantValue(self._param.Case.PrandtlNumber)


    def calculateSpeciesEnthalpy(self):
        """
        Compute and store species sensible enthalpies.

        This method constructs JANAF thermodynamic species for several common
        gas components and evaluates their sensible enthalpy expressions at
        the current temperature field ``T``. The resulting enthalpy fields are
        projected onto the P2 space and stored in the ``_hSpec`` dictionary.

        Notes
        -----
        The following species are currently supported:

        - CH₄
        - O₂
        - CO
        - CO₂
        - H₂O (gaseous)

        The resulting dictionary entries are:

        - ``_hSpec['CH4']``
        - ``_hSpec['O2']``
        - ``_hSpec['CO']``
        - ``_hSpec['CO2']``
        - ``_hSpec['H2O']``

        Each entry holds a finite element function representing the sensible
        enthalpy of the corresponding species.
        """
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
        Return vertex-based values of all mean-flow fields.

        The method converts the internal dictionary of finite element fields
        to vertex-based arrays and wraps them into a
        :class:`meanFlowVertexValues` instance.

        Returns
        -------
        meanFlowVertexValues
            Object containing vertex values for all fields in ``_fieldDict``,
            based on the mesh exported by ``self._mesh.exportMesh``.
        """
        return meanFlowVertexValues(self._fieldDict, self._mesh.exportMesh)

    def _getMeanFieldsToBeRead(self):
        """
        Assemble the list of mean-flow field names to be read from file.

        This method combines field names from several sources:

        - basic mean-flow fields specified in ``param``,
        - additional fields required by the energy handler,
        - additional fields required by the equation-of-state handler,
        - additional fields required by the reaction handler,
        - custom mean-flow quantities registered via
          :meth:`addCustomMeanFlowQuantity`.

        Duplicate names are removed to ensure that each field is read only
        once.

        Returns
        -------
        list of str
            Deduplicated list of field names to import from the mean-flow
            input file.

        Notes
        -----
        TODO: The tensor order is not currently encoded in the field names. This
        may be improved in the future to avoid hard-coded assumptions.
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
        """
        Register an additional mean-flow quantity to be read.

        Parameters
        ----------
        key : str
            Name of the additional mean-flow quantity that should be included
            when assembling the list of fields to import.

        Notes
        -----
        The key is appended to the internal list
        ``_customMeanFlowQuantities`` and will be included by
        :meth:`_getMeanFieldsToBeRead`. Duplicates are removed when the final
        list of field names is constructed.
        """
        self._customMeanFlowQuantities.append(key)

class meanFlowVertexValues(fieldProperties):
    """
    Vertex-based representation of mean-flow fields.

    This class converts `Field` objects storing finite element solutions
    into arrays of vertex values, providing a convenient interface for
    post-processing and exporting mean-flow quantities on the mesh vertices.

    **Initialize the meanFlowVertexValues object**

    Parameters
    ----------
    fieldDict : dict of {str: Field}
        Dictionary of mean-flow fields whose vertex values will be extracted.
    mesh : object
        Mesh object providing coordinates for the vertex arrays (typically
        ``mesh.exportMesh`` from the mean-flow class).

    Attributes
    ----------
    _fieldDict : dict
        Dictionary mapping field names to vertex-value arrays. For scalar
        fields, the values are stored as 1D arrays. For vector-valued fields,
        the values are stored as 2D arrays with shape
        ``(num_components, num_vertices)``.
    _isMean : bool
        Flag indicating that these values correspond to mean-flow quantities.
    _zeroField : ndarray
        Convenience array of zeros with length equal to the number of mesh
        vertices.
    _oneField : ndarray
        Convenience array of ones with length equal to the number of mesh
        vertices.

    Notes
    -----
    Vector-valued fields are handled by extracting each component from the
    corresponding subspace and assembling them into a single stacked array.
    """

    def __init__(self,fieldDict,mesh):
        """
        Initialize the meanFlowVertexValues instance.

        Parameters
        ----------
        fieldDict : dict of {str: Field}
            Dictionary containing the original finite element fields.
        mesh : object
            Mesh used to compute vertex coordinates and thus the vertex arrays.

        Notes
        -----
        For vector-valued fields (with more than one subspace), the method
        extracts coefficients for each component and stacks them into a
        two-dimensional array of shape
        ``(num_components, num_vertices)``. Scalar fields are stored as
        one-dimensional arrays of length ``num_vertices``.
        """
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
