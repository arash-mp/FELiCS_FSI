import  os
import  numpy               as np
from    .Field              import Field
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

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

    def __init__(self, FEMSpace, mesh, name=None, isStateVector=True, m=None, analysis=None):
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
        
        # Set defaults
        if not name:
            self.name = "q_hat"
        if not analysis:
            self.analysis = "Modal"
            logger.warning('No analysis type provided when initializing Mode object. Setting to "Modal" by default.')
    
        self.namesOfSubFields = self.getNamesOfSubFields()    
        
        self.isAdjoint  = False
        self.isResponse = False

    def getName(self):
        if self.isStateVector:
            return "q_hat"
        else:
            return super().getName()

    def setGain(self,gain):
        """
        Set the gain value for the mode.

        Parameters
        ----------
        gain : float
            Gain of the mode.
        """
        self._gain = gain

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


    def getGain(self):
        """
        Get the gain value of the mode.

        Returns
        -------
        float
            Gain of the mode. Returns -9999. if undefined.
        """
        try:
            return self._gain
        except: 
            logger.error('For this mode object no gain was defined. Returning "-9999."...')
            return -9999.


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
        ):
        
        # Get info on mode used to define filename
        if self.analysis == 'Modal':
            eigval      = self.getEigenValue()
            modeType    = 'Direct' if not self.isAdjoint else 'Adjoint'
            fileName    = f'{self.analysis}Solution_Omega_{modeType}_{np.round(eigval, 3)}.h5'
        elif self.analysis == 'Resolvent':
            frequency   = self.getFrequency()   
            modeType    = 'Response' if self.isResponse else 'Forcing'
            fileName    = f'{self.analysis}_Omega{np.round(frequency, 3)}_{modeType}_gain0.h5'   # NOTE: gain0 is a placeholder
        elif self.analysis == 'Input-Output':
            frequency   = self.getFrequency()
            modeType    = 'Response'
            fileName    = f'{self.analysis}_Omega{np.round(frequency, 3)}_{modeType}_gain0.h5'   # NOTE: Always gain 0 for IO modes
        else:
            logger.error(f'Analysis type "{self.analysis}" not recognized. Cannot import mode data.')
            return self, None
        
        # File name and group name
        importFilePath  = os.path.join(importDirPath, fileName)
        groupName       = "fluctuation/0/pointData/"  # TODO: remove all group names in FELiCS files
        
        # Call the reader from Field parent class
        self, notInFile = super().importData(
            reader,
            importFilePath,
            groupName
        )
        
        # Set the mode properties from the file
        
        
        return self, notInFile