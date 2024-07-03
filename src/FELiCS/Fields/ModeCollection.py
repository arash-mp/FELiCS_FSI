import numpy as np

from .Mode import Mode

from .fluctuationClass import fluctuationSolutions

class ModeCollection():

    def __init__(self, femSpace, mesh):
        self.modeList = []

        self.femSpace = femSpace
        self.mesh     = mesh


    def appendMode(self,mode):
        self.modeList.append(mode)


    def appendSolutionOfEigenProblem(self, solution, guess, adjoint=False):
        [eigVals, eigVecs, error] = solution
        numberOfModes = len(eigVals)

        for i in range(numberOfModes):
            mode = Mode(self.femSpace, self.mesh)
            mode.isAdjoint = adjoint
            mode.setError(error[i])
            mode.setGuess(guess)
            mode.setEigenValue(eigVals[i])
            mode.function.x.array[:] = eigVecs[i,:]
            self.modeList.append(mode)

    def getMaximumError(self):
        import numpy as np
        error = []
        for mode in self.modeList:
            error.append(mode.getError())
        if len(self.modeList)==0:
            pass #TODO: throw error
        else:
            return np.amax(error)

    def getDirectEigenValueSpectrum(self):
        spectrum = []
        for mode in self.modeList:
            if not mode.isAdjoint:
                spectrum.append(mode.getEigenValue())
        return spectrum


    def getNearestMode(self,gain=None,eigenValue=None,guess=None,waveNumber=None,frequency=None):
        #ToDo: get nearest mode to one of the above. Change handling of parameters
        pass


    def getLeadingMode(self,adjoint=False):
        leadingMode   = None
        if not adjoint:
            growthRateMax = -9990.
            for mode in self.modeList:
                if mode.isAdjoint == adjoint:
                    eigenValue = mode.getEigenValue()
                    if np.imag(eigenValue) > growthRateMax:
                        growthRateMax = np.imag(eigenValue)
                        leadingMode = mode
        elif adjoint:
            growthRateMin = 9990.
            for mode in self.modeList:
                if mode.isAdjoint == adjoint:
                    eigenValue = mode.getEigenValue()
                    if np.imag(eigenValue) < growthRateMin:
                        growthRateMin = np.imag(eigenValue)
                        leadingMode = mode
                    
        return leadingMode


    def getOldSolutionObject(self, meanFlow, param, FEMSpaces):
        # this is a wrapper for the old solution class, should be removed at the end of restructuring
        fluctSolutObjList    = []
        fluctSolutObjListAdj = []

        for mode in self.modeList:
            # construct for each EVal and EVec a fluctuationSolution
            fluctSolutObjList.append(fluctuationSolutions(
                                    param,
                                    meanFlow,
                                    FEMSpaces,
                                    mode.getEigenValue(),
                                    mode.function.x.array[:],
                                    mode.isAdjoint==False,
                                     )
            )

        return fluctSolutObjList


    def getSize(self):
        return len(self.modeList)

    def popList(self):
        return self.modeList.pop()
