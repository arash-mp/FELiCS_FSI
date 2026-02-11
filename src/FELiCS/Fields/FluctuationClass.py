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
# Standard libraries
from os import (
    listdir,
    remove,
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
from    FELiCS.Fields.FieldProperties                               import FieldProperties
from    FELiCS.Fields.Field                                         import Field
from    FELiCS.Equation.dependentVariables.EnergyHandler            import EnergyHandler
from    FELiCS.Equation.dependentVariables.EquationOfStateHandler   import EquationOfStateHandler
from    FELiCS.Equation.dependentVariables.HeatReleaseHandler       import HeatReleaseHandler
from    FELiCS.Equation.dependentVariables.MomentumHandler          import MomentumHandler
from    FELiCS.Equation.dependentVariables.ReactionHandler          import ReactionHandler
from    FELiCS.Misc.tensorUtils                                     import Tensor
from 	FELiCS.Misc.logging                                         import Logger

# Get the logger
logger = Logger.get_logger("felics")

class FluctuationClass(
    FieldProperties,
    ReactionHandler,
    EquationOfStateHandler,
    EnergyHandler,
    MomentumHandler,
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
    coordinate_system,
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

        # Useful objects
        self._param                 = param
        self._FEMSpaces             = FEMSpaces
        self._coordinateSystem      = coordinate_system
        FieldProperties.__init__(self)
        
        # Default flags
        self._isFluctuation         = True
        self._isMean                = False
        self._isSolution            = False
        
        # Zero-value tensors
        self._zeroVectorField       = Field(
        FEMSpaces.FunctionSpaceVectorVelocity,
        self._param.get_mesh(),
        name="zeroVector",
        ).get_tensor()
        self._zeroField             = Function(FEMSpaces.P2)
        self._zeroField             = Field(
        FEMSpaces.P2,
        self._param.get_mesh(),
        name="zeroScalar",
        ).get_tensor()
        self._fieldDict             = {}
        self._mean                  = mean
        self._transportedQuantities = param.get_transported_quantity_list()

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
            MomentumHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_momentum()

        if param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
            EquationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_eo_s()

        if param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
            EnergyHandler.__init__(
                                  self,
                                  )
            neededVariables += self._get_needed_fields_for_linear_energy()
            

        # Delete duplicates
        neededVariables = list(dict.fromkeys(neededVariables))

        # While not all needed fluctuations are calculated, try calculating them
        n_try = 1
        while not set(neededVariables).issubset((self._fieldDict.keys())):
            if param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
                self._relate_conservative_to_primitive_variables_momentum()
                self._initialize_molecular_momentum_diffusion_fluctuation()
            if param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
                self._initialize_eo_s_fluctuations()
            if param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
                self._relate_conservative_to_primitive_variables_energy()
                self._initialize_molecular_heat_diffusion_fluctuation()
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                logger.error('Attempt to calculate secondary variables not successful. Missing quantities: '\
                     + str(notInitializedFields) + ". Maybe the mixture file is still in the old format (ending with a '.mix' instead of '.json')?")

        # print(self._param.Mixture.getReactionMechanism()['type'])
        if not self._param.Mixture.get_reaction_mechanism()['type'] == 'None':
            ReactionHandler.__init__(
                self,
                )
            self._initialize_reactions()
            
class FluctuationSolutions(
    FieldProperties,
    HeatReleaseHandler,
    EquationOfStateHandler,
    MomentumHandler,
    EnergyHandler,
    ReactionHandler,
):
    """
    Stores and manages linearized fluctuation solutions.

    This class contains the computed fluctuation solutions in the mixed
    function space. 

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
    gain_number=-1,
    gain_value=-1,
    
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
        self._gainValue = np.real(gain_value)
        self._omega = omega
        self._isResponseOrDirect = isResponseOrDirect
        self._vmixedVector = vmixedVector
        self._gainNumber = gain_number
        self._FEMSpaces = FEMSpaces
        self._mean = mean
        # self._fieldDict = mean.fieldDict
        self._transportedQuantities = param.get_transported_quantity_list()
        self._param = param

    def _fluc_export_wrapper(self,
    group,
    ):
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
            MomentumHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_momentum()

        if self._param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
            EquationOfStateHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_eo_s()

        if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
            EnergyHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_energy()
    
        if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
            EnergyHandler.__init__(
                self,
                )
            neededVariables += self._get_needed_fields_for_linear_energy()
   

        # Delete duplicates
        neededVariables = list(dict.fromkeys(neededVariables))

        # While not all needed fluctuations are calculated, try calculating them
        n_try = 0
        while not set(neededVariables).issubset((self._fieldDict.keys())):
            meanVertexValues = self._mean.get_vertex_values()
            if self._param.Case.SetOfEquations['Momentum']['Equation'] not in ['None']:
                self._relate_conservative_to_primitive_variables_momentum(
                                                    meanVertexValues,
                                                    )
                self._initialize_molecular_momentum_diffusion_fluctuation(meanVertexValues)

            if self._param.Case.SetOfEquations['EquationOfState']['Equation'] not in ['None']:
                self._initialize_eo_s_fluctuations(
                                            meanVertexValues,
                                                )

            if self._param.Case.SetOfEquations['Energy']['Equation'] not in ['None']:
                self._relate_conservative_to_primitive_variables_energy(
                                            meanVertexValues,
                                                )
                self._initialize_molecular_heat_diffusion_fluctuation(
                                            meanVertexValues,
                                            )
            n_try += 1
            if n_try > 100:
                notInitializedFields = list(set(neededVariables) - set(list(self._fieldDict.keys())))
                logger.error('Attempt to calculate secondary variables not successful. Missing quantities: ' +notInitializedFields +". Maybe the mixture file is still in the old format (ending with a '.mix' instead of '.json')?")

        if not self._param.Mixture.get_reaction_mechanism()['type'] == 'None':
            ReactionHandler.__init__(
                self,
                )
            self._initialize_reactions(
                meanVertexValues,
                )
        self._writeDictToH5(
        self._fieldDict,
        group,
        True,
        )

    def _import_sol_vector(self,
    filename,
    ):
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
 
        hf = File(
        filename,
        'r',
        )
        frequency = hf['fluctuation/0/'].attrs.get('frequency')
        fieldMagnitude = np.array(hf[f'fluctuation/0/{frequency}/magnitude'][:])
        fieldAngle = np.array(hf[f'fluctuation/0/{frequency}/angle'][:])

        return fieldMagnitude * np.exp(1j * fieldAngle)


    @property
    def solut_vector(self):
        """
        Complex-valued vector of the fluctuation solution.

        Returns
        -------
        np.ndarray
            Solution vector in mixed function space.
        """
        return self._vmixedVector


    @property
    def solution_kind(self):
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
    def gain_number(self):
        """
        Gain number identifier.

        Returns
        -------
        int
            Gain index or -1 if unused.
        """
        return self._gainNumber


    @property
    def gain_value(self):
        """
        Gain magnitude.

        Returns
        -------
        float
            Real-valued gain.
        """
        return self._gainValue

