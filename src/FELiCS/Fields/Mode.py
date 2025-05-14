
from 	FELiCS.Misc.logging import Logger
from    .Field              import Field

# Get the logger
logger = Logger.get_logger("felics")

class Mode(Field):

    def __init__(self,FEMSpace, mesh):
        super().__init__(FEMSpace, mesh)
        self.isAdjoint  = False
        self.isResponse = False

    def setGain(self,gain):
        self._gain = gain

    def setFrequency(self,frequency):
        self._frequency = frequency

    def setEigenValue(self,eigenValue):
        self._eigenValue = eigenValue


    def setWaveNumber(self,waveNumber):
        self._waveNumber = waveNumber


    def setGuess(self,guess):
        self._guess = guess
        
    def setError(self,error):
        self._error = error


    def getGain(self):
        try:
            return self._gain
        except: 
            logger.error('For this mode object no gain was defined. Returning "-9999."...')
            return -9999.


    def getFrequency(self):
        try:
            return self._frequency
        except: 
            logger.error('For this mode object no frequency was defined. Returning "-9999."...')
            return -9999.


    def getEigenValue(self):
        try:
            return self._eigenValue
        except: 
            logger.error('For this mode object no eigen value was defined. Returning "-9999."...')
            return -9999.


    def getWaveNumber(self):
        try:
            return self._waveNumber
        except: 
            logger.error('For this mode object no waveNumber was defined. Returning "-9999."...')
            return -9999.


    def getGuess(self):
        try:
            return self._guess
        except: 
            logger.error('For this mode object no guess was defined. Returning "-9999."...')
            return -9999.

    def getError(self):
        try:
            return self._error
        except: 
            logger.error('For this mode object no error was defined. Returning "-9999."...')
            return -9999.


