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
            mode.setError(error)
            mode.setGuess(guess)
            mode.setEigenValue(eigVals[i])
            mode.function.x.array[:] = eigVecs[i,:]
            self.modeList.append(mode)

    def getMaximumError(self):
        import numpy as np
        error = []
        for mode in self.modeList:
            error.append(mode.getError())
        return np.amax(error)


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
