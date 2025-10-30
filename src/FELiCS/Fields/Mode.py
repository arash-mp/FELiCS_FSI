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

    Attributes
    ----------
    isAdjoint : bool
        Indicates if the mode is an adjoint mode.
    isResponse : bool
        Indicates if the mode is a response mode.
    """

    def __init__(self, FEMSpace, mesh, name="q_hat", isStateVector=True, m=0, analysisType='Modal'):
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
        # NOTE: why are these hardcoded here?
        
        # Set values
        self.analysisType      = AnalysisType[analysisType.upper()]
    
        # Define name of the subfields (variables of the mode)
        self.namesOfSubFields   = self.getNamesOfSubFields()

        # Set some defaults in not a good way .> TODO: fix this as a property
        if self.analysisType == AnalysisType.MODAL:
            self.modeType = ModeType.DIRECT
        elif self.analysisType == AnalysisType.RESOLVENT:
            self.modeType = ModeType.FORCING
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            self.modeType = ModeType.RESPONSE
         
        self.isAdjoint          = False
        self.isResponse         = False

    # TODO: fix the setter/getter methods with properties

    @property
    def name(self):
        if self.isStateVector:
            return "q_hat"
        else:
            return super().name


    @property
    def gain(self):
        try:
            return self._gain
        except: 
            logger.error('For this mode object no gain was defined. Returning "-9999."...')
            return -9999.

    @gain.setter
    def gain(self, gain):
        self._gain = gain



    def setGainNumber(self,gainNumber):
        """
        Set the gain number for the mode.

        Parameters
        ----------
        gainNumber : int
            Gain number of the mode.
        """
        self._gainNumber = gainNumber

    def setFrequency(self,frequency):
        """
        Set the frequency for the mode.

        Parameters
        ----------
        frequency : float
            Frequency associated with the mode.
        """
        self._frequency = frequency

    def setEigenValue(self,eigenValue):
        """
        Set the eigenvalue for the mode.

        Parameters
        ----------
        eigenValue : float
            Eigenvalue associated with the mode.
        """
        self._eigenValue = eigenValue


    def setWaveNumber(self,waveNumber):
        """
        Set the wave number for the mode.
        # NOTE: Why do we have this AND self.m?

        Parameters
        ----------
        waveNumber : float
            Wave number corresponding to the mode.
        """
        self._waveNumber = waveNumber


    def setGuess(self,guess):
        """
        Set the initial guess for the mode.

        Parameters
        ----------
        guess : float
            Initial guess used in the mode computation.
        """
        self._guess = guess
        
    def setError(self,error):
        """
        Set the error value for the mode.

        Parameters
        ----------
        error : float
            Error associated with the mode solution.
        """
        self._error = error


    def getGainNumber(self):
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
            logger.error('For this mode object no gain number was defined. Returning "-1"...')
            return -1


    def getFrequency(self):
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
            logger.error('For this mode object no frequency was defined. Returning "-9999."...')
            return -9999.


    def getEigenValue(self):
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
            logger.error('For this mode object no eigen value was defined. Returning "-9999."...')
            return -9999.


    def getWaveNumber(self):
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
            logger.error('For this mode object no waveNumber was defined. Returning "-9999."...')
            return -9999.


    def getGuess(self):
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
            logger.error('For this mode object no guess was defined. Returning "-9999."...')
            return -9999.

    def getError(self):
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
            logger.error('For this mode object no error was defined. Returning "-9999."...')
            return -9999.
        
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
                eigval      = self.getEigenValue()
                modeType    = 'Direct' if not self.isAdjoint else 'Adjoint'
                fileName    = f'ModalSolution_Omega_{modeType}_{np.round(eigval, 3)}.h5'
            elif self.analysisType == AnalysisType.RESOLVENT:
                frequency   = self.getFrequency()   
                modeType    = 'Response' if self.isResponse else 'Forcing'
                fileName    = f'Resolvent_Omega{np.round(frequency, 3)}_{modeType}_gain{self.getGainNumber()}.h5'
            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                frequency   = self.getFrequency()
                fileName    = f'Input-Output_Omega{np.round(frequency, 3)}_Response_gain0.h5'   # NOTE: Always gain 0 for IO modes
        
        # File name and group name
        importFilePath      = os.path.join(importDirPath, fileName)
        groupName           = "fluctuation/0/pointData/"  # TODO: remove all group names in FELiCS files
        
        # Call the reader from Field parent class
        self, notInFile     = super().importData(
            reader,
            importFilePath,
            groupName
        )
        
        # Read eignvalue or gain from file
        if self.analysisType == AnalysisType.MODAL:
            with h5py.File(importFilePath, 'r') as f:
                eigval = complex(f["fluctuation/0"].attrs['frequency']) # TODO: save the eigenvalue not as a string in files!
                self.setEigenValue(eigval)
        
        elif self.analysisType in [AnalysisType.RESOLVENT, AnalysisType.INPUT_OUTPUT]:
            with h5py.File(importFilePath, 'r') as f:
                freq_string = f["fluctuation/0"].attrs['frequency']
                if 'j' in freq_string:
                    frequency = complex(freq_string)
                else:
                    frequency = float(freq_string)
                self.setFrequency(frequency)
        
        # TODO: Save the gain and eigenvalues in files and set them here
        
        return self, notInFile
