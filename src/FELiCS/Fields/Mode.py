import  os
import  h5py
import  numpy               as np
from    enum                import Enum
from    .Field              import Field
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class ModeType(Enum):
    NONE     = 0
    DIRECT   = 1
    ADJOINT  = 2
    RESPONSE = 3
    FORCING  = 4

class AnalysisType(Enum):
    NONE         = 0
    MODAL        = 1
    RESOLVENT    = 2
    INPUT_OUTPUT = 3


class Mode(Field):
    """
    Field-based representation of a computational mode.

    This class extends the Field object to represent modes in a FEM-based
    analysis, such as eigenmodes or response modes. It includes properties
    like gain, frequency, eigenvalue, wave number, and error metrics.

    **Initialize the Mode object**

    Parameters
    ----------
    FEMSpace : object
        The finite element space defining the discretization.
    mesh : object
        The mesh on which the FEM space is defined.

    """

    def __init__(self, FEMSpace, mesh, name="q_hat", isStateVector=True, m=0, analysisType='Modal', modeType = None):
        """
        Initializes the Mode instance.

        Parameters
        ----------
        FEMSpace : object
            The finite element space defining the discretization.
        mesh : object
            The mesh on which the FEM space is defined.
        """
        super().__init__(FEMSpace, mesh, name, isStateVector, m)
        
        # Set values
        self.analysisType      = AnalysisType[analysisType.upper()]
        # TODO Sophie: write error if analysisType is not given
    
        # Define name of the subfields (variables of the mode)
        self.namesOfSubFields   = self.getNamesOfSubFields()

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
    def name(self):
        if self.isStateVector:
            return "q_hat"
        else:
            return super().name


    @property
    def omega(self):
        if self.analysisType == AnalysisType.MODAL:
            return self.eigenValue
        elif self.analysisType in [AnalysisType.RESOLVENT, AnalysisType.INPUT_OUTPUT] :
            return self.frequency

    @property
    def gain(self):
        try:
            return self._gain
        except: 
            logger.warning('For this mode object no gain was defined. Returning "-9999."...')
            return -9999.

    @gain.setter
    def gain(self, gain):
        self._gain = gain



    @property
    def gainNumber(self):
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

    @gainNumber.setter
    def gainNumber(self, gainNumber):
        """
        Set the gain number for the mode.

        Parameters
        ----------
        gainNumber : int
            Gain number of the mode.
        """
        self._gainNumber = gainNumber

    @property
    def frequency(self):
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
    def frequency(self,frequency):
        """
        Set the frequency for the mode.

        Parameters
        ----------
        frequency : float
            Frequency associated with the mode.
        """
        self._frequency = frequency


    @property
    def eigenValue(self):
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

    @eigenValue.setter
    def eigenValue(self,eigenValue):
        """
        Set the eigenvalue for the mode.

        Parameters
        ----------
        eigenValue : float
            Eigenvalue associated with the mode.
        """
        self._eigenValue = eigenValue

    @property
    def waveNumber(self):
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
 
    @waveNumber.setter
    def waveNumber(self,waveNumber):
        """
        Set the wave number for the mode.
        # NOTE: Why do we have this AND self.m?

        Parameters
        ----------
        waveNumber : float
            Wave number corresponding to the mode.
        """
        self._waveNumber = waveNumber

    @property
    def guess(self):
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
    def guess(self,guess):
        """
        Set the initial guess for the mode.

        Parameters
        ----------
        guess : float
            Initial guess used in the mode computation.
        """
        self._guess = guess

    @property 
    def error(self):
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
    def error(self,error):
        """
        Set the error value for the mode.

        Parameters
        ----------
        error : float
            Error associated with the mode solution.
        """
        self._error = error
        
    def describe(self):
        """
        Print a description of the mode, including its properties.
        """
        logger.info(f"  Mode Type:     {self.modeType.name}")
        if self.analysisType == AnalysisType.MODAL:
            logger.info(f"  Guess:         {self.guess}")
            logger.info(f"  Eigenvalue:    {self.eigenValue}")
        elif self.analysisType == AnalysisType.RESOLVENT:
            logger.info(f"  Frequency:     {self.frequency}")
            logger.info(f"  Gain:          {self.gain}")
            logger.info(f"  Gain Number:   {self.gainNumber}")
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            logger.info(f"  Frequency:     {self.frequency}")
        logger.info(f"  Wave Number:   {self.waveNumber}") 

    def exportToH5(self, writer, fileName=None):
        # create standard fileName if none is given
        if fileName is None:
            fileName = "Mode_" \
                       + self.analysisType.name.capitalize() + "_" \
                       + self.modeType.name.capitalize() + "_" \
                       + "Omega_" \
                       + "{:.3f}".format(self.omega)
            if self.analysisType == AnalysisType.RESOLVENT:
                fileName += "_GainNb_" + str(self.gainNumber)

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
            attrGainNumber.value = self.gainNumber
            attrList.append(attrGainNumber)

        # export mode        
        writer.exportFieldToH5(self, fileName, attrList)


    def importData(
            self,
            reader,
            importDirPath,
            importFileName = None,
        ):
        
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
                fileName    = f'Mode_Resolvent_{modeType}_Omega_{"{:.3f}".format(self.frequency)}_GainNb_{self.gainNumber}.h5'
            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                fileName    = f'Mode_Input_output_Response_Omega_{"{:.3f}".format(self.frequency)}.h5'
        
        # File name
        importFilePath      = os.path.join(importDirPath, fileName)
        
        # Call the reader from Field parent class
        self, notInFile     = super().importData(
            reader,
            importFilePath,
        )
        
        # Read eignvalue or gain from file
        if self.analysisType == AnalysisType.MODAL:
            with h5py.File(importFilePath, 'r') as f:
                self.eigenValue = f["omega"][()]
                
        elif self.analysisType == AnalysisType.RESOLVENT:
            with h5py.File(importFilePath, 'r') as f:
                self.frequency  = f["omega"][()]     # NOTE: slight inconsistency in naming
                self.gain       = f["gain"][()]
                self.gainNumber = f["number"][()]
                
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            with h5py.File(importFilePath, 'r') as f:
                self.frequency  = f["omega"][()]    # NOTE: slight inconsistency in naming
        
        return self, notInFile
