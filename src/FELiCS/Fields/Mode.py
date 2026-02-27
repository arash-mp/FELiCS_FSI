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
from    enum import Enum
import  os

# Third party libraries
import h5py
import numpy as np

# Local Libraries and methods
from FELiCS.Fields.Field import Field
from FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class ModeType(Enum):
    """
    Enumeration of mode roles in the analysis.

    Members
    -------
    NONE : 0
        No specific mode type.
    DIRECT : 1
        Direct eigenmode.
    ADJOINT : 2
        Adjoint eigenmode.
    RESPONSE : 3
        Response mode in resolvent or input–output analysis.
    FORCING : 4
        Forcing mode in resolvent or input–output analysis.
    """
    NONE     = 0
    DIRECT   = 1
    ADJOINT  = 2
    RESPONSE = 3
    FORCING  = 4

class AnalysisType(Enum):
    """
    Enumeration of analysis types supported for modes.

    Members
    -------
    NONE : 0
        No analysis type.
    MODAL : 1
        Modal (eigenvalue) analysis.
    RESOLVENT : 2
        Resolvent analysis.
    INPUT_OUTPUT : 3
        Input–output analysis.
    """
    NONE         = 0
    MODAL        = 1
    RESOLVENT    = 2
    INPUT_OUTPUT = 3


class Mode(Field):
    """
    Field-based representation of a computational mode.

    This class extends :class:`Field` to represent modes in a FEM-based
    analysis, such as eigenmodes, forcing modes or response modes. It
    stores additional metadata such as gain, frequency, eigenvalue, wave
    number, and error metrics, together with an analysis type and a mode
    type.

    **Initialize the Mode object**

    Parameters
    ----------
    FEMSpace : dolfinx.fem.FunctionSpace
        The finite element space defining the discretization.
    mesh : FELiCS.SpaceDisc.FELiCSMesh
        Mesh on which the FEM space is defined.
    name : str, optional
        Name of the underlying field. Default is ``"q_hat"``.
    isStateVector : bool, optional
        If True, this mode is treated as a state vector in a mixed
        formulation. Default is True.
    m : int, optional
        Wave number associated with a spectral spatial dimension.
        Default is 0.
    analysisType : str, optional
        Type of analysis for which this mode is defined. The string is mapped
        (case-insensitively) to :class:`AnalysisType` (e.g. ``'Modal'``,
        ``'Resolvent'``). Default is ``'Modal'``.
    modeType : str or None, optional
        Mode role within the chosen analysis type, mapped to
        :class:`ModeType` (e.g. ``'direct'``, ``'adjoint'``, ``'response'``,
        ``'forcing'``). If ``None``, a sensible default is chosen based on
        ``analysisType``.

    Attributes
    ----------
    analysisType : AnalysisType
        The specific analysis type enum associated with this mode.
    modeType : ModeType
        The specific role of the mode (e.g. DIRECT, ADJOINT, FORCING).
    waveNumber : int
        Wave number associated with the mode (initialized by `m`).
    namesOfSubFields : list
        List of names of the subfields (variables) constituting the mode.
    omega : float
        Property returning the eigenvalue (Modal) or frequency (Resolvent/IO).
    gain : float
        Property returning the mode gain (returns -9999. if undefined).
    error : float
        Property returning the error associated with the mode (returns -9999. if undefined).
    """

    def __init__(
        self,
        FEMSpace,
        mesh,
        name="q_hat",
        isStateVector=True,
        m=0,
        analysisType='Modal',
        modeType = None
    ):
        """
        Initialize a Mode instance.

        Parameters
        ----------
        FEMSpace : dolfinx.fem.FunctionSpace
            The finite element space defining the discretization.
        mesh : FELiCS.SpaceDisc.FELiCSMesh
            Mesh on which the FEM space is defined.
        name : str, optional
            Name of the underlying field. Default is ``"q_hat"``.
        isStateVector : bool, optional
            If True, this mode is treated as a state vector in a mixed
            formulation. Default is True.
        m : int, optional
            Wave number associated with a spectral spatial dimension.
            Default is 0.
        analysisType : {'Modal', 'Resolvent', 'Input_Output'}, optional
            Analysis type. The string is converted to an
            :class:`AnalysisType` enum. Default is ``'Modal'``.
        modeType : str or None, optional
            Mode role within the chosen analysis type, mapped to
            :class:`ModeType`. If ``None``, a default is selected based on
            ``analysisType`` (e.g. DIRECT for MODAL, FORCING for RESOLVENT,
            RESPONSE for INPUT_OUTPUT).
        """
        super().__init__(
            FEMSpace,
            mesh,
            name,
            isStateVector,
            m
        )
        
        # Set values
        self.analysisType      = AnalysisType[analysisType.upper()]
        # TODO Sophie: write error if analysisType is not given
        self.wave_number        = m
    
        # Define name of the subfields (variables of the mode)
        self.namesOfSubFields   = self.get_names_of_sub_fields()

        if modeType is None:
            # Set some defaults in not a good way .> TODO: fix this as a property
            if self.analysisType == AnalysisType.MODAL:
                self.modeType = ModeType.DIRECT
            elif self.analysisType == AnalysisType.RESOLVENT:
                self.modeType = ModeType.FORCING
            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                self.modeType = ModeType.RESPONSE
        else:
            self.modeType = ModeType[modeType.upper()]


    @property
    def name(
        self
    ):
        """
        Name of the mode.

        For state-vector modes (``isStateVector`` is True), this is always
        returned as ``"q_hat"``, irrespective of the underlying field name.
        Otherwise, the name from the base :class:`Field` class is used.
        """
        if self.isStateVector:
            return "q_hat"
        else:
            return super().name


    @property
    def omega(
        self
    ):
        """
        Angular frequency associated with the mode.

        For modal analysis, this is the eigenvalue. For resolvent and
        input–output analysis, this is the real frequency.

        Returns
        -------
        float
            Eigenvalue (modal) or frequency (resolvent / input–output).
        """
        if self.analysisType == AnalysisType.MODAL:
            return self.eigen_value
        elif self.analysisType in [AnalysisType.RESOLVENT, AnalysisType.INPUT_OUTPUT] :
            return self.frequency

    @property
    def gain(
        self
    ):
        """
        Gain associated with the mode.

        Returns
        -------
        float
            The mode gain. If undefined, returns ``-9999.`` and issues a warning.
        """
        try:
            return self._gain
        except: 
            logger.warning('For this mode object no gain was defined. Returning "-9999."...')
            return -9999.

    @gain.setter
    def gain(
        self,
        gain
    ):
        """
        Set the gain of the mode.

        Parameters
        ----------
        gain : float
            Gain value associated with the mode.
        """
        self._gain = gain



    @property
    def gain_number(
        self
    ):
        """
        Get the gain number of the mode.

        Returns
        -------
        int
            Gain number of the mode. Returns -1 if undefined.
        """
        try:
            return self._gainNumber
        except: 
            logger.warning('For this mode object no gain number was defined. Returning "-1"...')
            return -1

    @gain_number.setter
    def gain_number(
        self,
        gain_number
    ):
        """
        Set the gain number for the mode.

        Parameters
        ----------
        gainNumber : int
            Gain number of the mode.
        """
        self._gainNumber = gain_number

    @property
    def frequency(
        self
    ):
        """
        Get the frequency of the mode.

        Returns
        -------
        float
            Frequency of the mode. Returns -9999. if undefined.
        """
        try:
            return self._frequency
        except: 
            logger.warning('For this mode object no frequency was defined. Returning "-9999."...')
            return -9999.

    @frequency.setter
    def frequency(
        self,
        frequency
    ):
        """
        Set the frequency for the mode.

        Parameters
        ----------
        frequency : float
            Frequency associated with the mode.
        """
        self._frequency = frequency


    @property
    def eigen_value(
        self
    ):
        """
        Get the eigenvalue of the mode.

        Returns
        -------
        float
            Eigenvalue of the mode. Returns -9999. if undefined.
        """
        try:
            return self._eigenValue
        except: 
            logger.warning('For this mode object no eigen value was defined. Returning "-9999."...')
            return -9999.

    @eigen_value.setter
    def eigen_value(
        self,
        eigen_value
    ):
        """
        Set the eigenvalue for the mode.

        Parameters
        ----------
        eigenValue : float
            Eigenvalue associated with the mode.
        """
        self._eigenValue = eigen_value

    @property
    def wave_number(
        self
    ):
        """
        Get the wave number of the mode.

        Returns
        -------
        float
            Wave number of the mode. Returns -9999. if undefined.
        """
        try:
            return self._waveNumber
        except: 
            logger.warning('For this mode object no waveNumber was defined. Returning "-9999."...')
            return -9999.
 
    @wave_number.setter
    def wave_number(
        self,
        wave_number
    ):
        """
        Set the wave number for the mode.

        Parameters
        ----------
        waveNumber : float or int
            Wave number corresponding to the mode.
        """
        self._waveNumber = wave_number

    @property
    def guess(
        self
    ):
        """
        Get the initial guess of the mode.

        Returns
        -------
        float
            Initial guess used. Returns -9999. if undefined.
        """
        try:
            return self._guess
        except: 
            logger.warning('For this mode object no guess was defined. Returning "-9999."...')
            return -9999.

    @guess.setter
    def guess(
        self,
        guess
    ):
        """
        Set the initial guess for the mode.

        Parameters
        ----------
        guess : float
            Initial guess used in the mode computation.
        """
        self._guess = guess

    @property 
    def error(
        self
    ):
        """
        Get the error associated with the mode.

        Returns
        -------
        float
            Error of the mode. Returns -9999. if undefined.
        """
        try:
            return self._error
        except: 
            logger.warning('For this mode object no error was defined. Returning "-9999."...')
            return -9999.

    @error.setter
    def error(
        self,
        error
    ):
        """
        Set the error value for the mode.

        Parameters
        ----------
        error : float
            Error associated with the mode solution.
        """
        self._error = error
        
    def describe(
        self
    ):
        """
        Log a summary of the mode's properties.

        Prints the mode type, and depending on the analysis type, logs
        relevant metrics such as guess and eigenvalue (Modal) or frequency,
        gain, and gain number (Resolvent/Input-Output).
        """
        logger.info(f"  Mode Type:     {self.modeType.name}")
        if self.analysisType == AnalysisType.MODAL:
            logger.info(f"  Guess:         {self.guess}")
            logger.info(f"  Eigenvalue:    {self.eigen_value}")
        elif self.analysisType == AnalysisType.RESOLVENT:
            logger.info(f"  Frequency:     {self.frequency}")
            logger.info(f"  Gain:          {self.gain}")
            logger.info(f"  Gain Number:   {self.gain_number}")
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            logger.info(f"  Frequency:     {self.frequency}")
        logger.info(f"  Wave Number:   {self.wave_number}") 

    def plot(
        self,
        variableName=None,
        xlim=None,
        ylim=None,
        plotType="real",
        clim=None
    ):
        """
        Plot a component of the mode using the base Field plot method.

        For mixed modes, a specific variable can be selected by name. If no
        name is provided, the first scalar field or the first component of
        the first vector field is plotted.

        Parameters
        ----------
        variableName : str or None, optional
            Name of the variable to plot. This may refer to a mixed subfield
            name or a component of a vector subfield (e.g. ``u_x``). If None,
            a default component is chosen.
        xlim, ylim : tuple, optional
            Axis limits passed to the underlying Field plotter.
        plotType : {'real', 'imag', 'magnitude'}, optional
            Component of the complex field to plot. Default is ``"real"``.
        clim : tuple, optional
            Color limits passed to the underlying Field plotter.
        """

        
        # Check that the mode is defined on a 2D mesh, 3D not yet implemented
        if self.mesh.dolfinxMesh.topology.dim != 2:
            logger.warning("Mode.plot(): plotting is currently only implemented for 2D meshes. Returning without plotting.")
            return

        # If the mode is already scalar, plot directly
        if self.info['type'] == 'scalar':
            return super().plot(
                                xlim=xlim,
                                ylim=ylim,
                                plotType=plotType,
                                clim=clim
                                )

        fields  = self.getListOfSubFields()
        names   = self.getNamesOfSubFields()

        selected_field = None

        if variableName is not None:
            if variableName in names:
                selected_field = fields[names.index(variableName)]
            else:

                # Try to match a vector component name inside vector subfields
                for field in fields:
                    if field.info['type'] == 'vector':
                        component_names     = field.getNamesOfSubFields()
                        if variableName in component_names:
                            selected_field  = field.getListOfSubFields()[
                                component_names.index(variableName)
                            ]
                            break
                if selected_field is None:
                    logger.warning(
                        f"Mode.plot(): variable '{variableName}' not found. "
                        "Using the default component instead."
                    )

        if selected_field is None:

            # Default: first scalar field or first component of first vector

            first_field = fields[0]
            if first_field.info['type']     == 'vector':
                selected_field = first_field.getListOfSubFields()[0]
            elif first_field.info['type']   == 'mixed':

                # Take the first scalar component from the mixed subfield
                nested_fields   = first_field.getListOfSubFields()
                nested_first    = nested_fields[0]
                if nested_first.info['type'] == 'vector':
                    selected_field  = nested_first.getListOfSubFields()[0]
                else:
                    selected_field  = nested_first
            else:
                selected_field      = first_field

        return selected_field.plot(
                                    xlim=xlim,
                                    ylim=ylim,
                                    plotType=plotType,
                                    clim=clim
                                   )

    def export_to_h5(
        self,
        writer,
        fileName=None
     ):

        """
        Export the mode to an HDF5 file.

        If no file name is provided, a default name is constructed from the
        analysis type, mode type, and frequency/eigenvalue (and gain number
        for resolvent modes). Mode metadata such as omega and, where
        applicable, gain and gain number are written as attributes.

        Parameters
        ----------
        writer : object
            Writer object providing an
            ``exportFieldToH5(field, fileName, attrList)`` method.
        fileName : str, optional
            Base file name (without directory). If ``None``, a default is
            generated based on analysis type and mode properties.
        """

        # create standard fileName if none is given
        if fileName is None:
            fileName = ("Mode_" 
                + self.analysisType.name.capitalize() + "_" 
                + self.modeType.name.capitalize() + "_" 
                + "Omega_" 
                + "{:.3f}".format(self.omega)
            )
            if self.analysisType == AnalysisType.RESOLVENT:
                fileName += "_GainNb_" + str(self.gain_number)

        # give a list of attributes to store in the file
        attr = lambda : None
        attr.name  = "omega"
        attr.value = self.omega
        attrList = [attr]
        if self.analysisType in [AnalysisType.RESOLVENT, AnalysisType.INPUT_OUTPUT]:
            attrGain  = lambda : None
            attrGain.name = "gain"
            attrGain.value = self.gain
            attrList.append(attrGain)

        if self.analysisType in [AnalysisType.RESOLVENT]:
            attrGainNumber  = lambda : None
            attrGainNumber.name = "number"
            attrGainNumber.value = self.gain_number
            attrList.append(attrGainNumber)

        # export mode        
        writer.export_field_to_h5(
            self,
            fileName,
            attrList,
        )


    def import_data(
        self,
        reader,
        importDirPath,
        importFileName = None
    ):
        """
        Import mode data from an HDF5 file.

        If no explicit file name is provided, a default name is constructed
        from the analysis type, mode type, and frequency/eigenvalue (and
        gain number for resolvent modes). After importing the field
        coefficients via the parent :class:`Field` method, the eigenvalue,
        frequency, gain and gain number are read back from file attributes.

        Parameters
        ----------
        reader : object
            Reader object providing an ``importInField(field, filePath, groupName)``
            method.
        importDirPath : str
            Directory from which the mode file is read.
        importFileName : str, optional
            Explicit file name to import from. If ``None``, a default name
            is constructed from the mode metadata.

        Returns
        -------
        Mode
            The updated mode (self) with imported coefficients and metadata.
        list of str
            List of variable names that were not found in the file.

        Notes
        -----
        For Resolvent and Input-Output analyses, the mode frequency is populated
        from the 'omega' dataset in the HDF5 file.
        """
        # In case we want to import from a specific file
        if importFileName is not None:
            fileName        = importFileName
        else:
            # Check that analysis type is set to set a default name
            if self.analysisType == AnalysisType.MODAL:
                modeType    = 'Direct' if self.modeType is ModeType.DIRECT else 'Adjoint'
                fileName    = f'Mode_Modal_{modeType}_Omega_{"{:.3f}".format(self.omega)}.h5'
            elif self.analysisType == AnalysisType.RESOLVENT:
                modeType    = 'Response' if self.modeType is ModeType.RESPONSE else 'Forcing'
                fileName    = f'Mode_Resolvent_{modeType}_Omega_{"{:.3f}".format(self.frequency)}_GainNb_{self.gain_number}.h5'
            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                fileName    = f'Mode_Input_output_Response_Omega_{"{:.3f}".format(self.frequency)}.h5'
        
        # File name
        importFilePath      = os.path.join(
            importDirPath,
            fileName,
        )
        
        # Call the reader from Field parent class
        self, notInFile     = super().import_data(
            reader,
            importFilePath,
        )
        
        # Read eignvalue or gain from file
        if self.analysisType == AnalysisType.MODAL:
            with h5py.File(
                importFilePath,
                'r',
            ) as f:
                self.eigen_value = f["omega"][()]
                
        elif self.analysisType == AnalysisType.RESOLVENT:
            with h5py.File(
                importFilePath,
                'r',
            ) as f:
                self.frequency  = f["omega"][()]     # NOTE: slight inconsistency in naming
                self.gain       = f["gain"][()]
                self.gain_number = f["number"][()]
                
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            with h5py.File(
                importFilePath,
                'r',
            ) as f:
                self.frequency  = f["omega"][()]    # NOTE: slight inconsistency in naming
        
        return self, notInFile
