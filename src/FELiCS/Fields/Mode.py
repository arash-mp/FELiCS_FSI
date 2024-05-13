
from FELiCS.Misc.functions import printError

from .Field import Field

class Mode(Field):

    def __init__(FEMSpace, mesh):
        super().__init(FEMSpace, mesh)
        # omega /lambda
        # guess
        #  


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
        


    def getGain(self):
        try:
            return self._gain
        except: 
            printError('For this mode object no gain was defined. Returning 0..')
            return 0.


    def getFrequency(self):
        try:
            return self._frequency
        except: 
            printError('For this mode object no frequency was defined. Returning 0..')
            return 0.


    def getEigenValue(self):
        try:
            return self._eigenValue
        except: 
            printError('For this mode object no eigen value was defined. Returning 0..')
            return 0.


    def getWaveNumber(self):
        try:
            return self._waveNumber
        except: 
            printError('For this mode object no waveNumber was defined. Returning 0..')
            return 0.


    def getGuess(self):
        try:
            return self._guess
        except: 
            printError('For this mode object no guess was defined. Returning 0..')
            return 0.


