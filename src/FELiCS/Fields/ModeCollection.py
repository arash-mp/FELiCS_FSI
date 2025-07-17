import numpy as np

from .Mode import Mode

from .fluctuationClass import fluctuationSolutions

class ModeCollection():
    """
    Container for managing a collection of Mode objects.

    This class provides methods to append, retrieve, and analyze modes, as well as
    to interface with legacy solution objects. It is designed to store and manage
    eigenmodes and their properties for further analysis.

    **Initialize the ModeCollection object**

    Parameters
    ----------
    femSpace : object
        The finite element space associated with the modes.
    mesh : object
        The mesh associated with the modes.

    Attributes
    ----------
    modeList : list
        List of Mode objects in the collection.
    femSpace : object
        The finite element space associated with the modes.
    mesh : object
        The mesh associated with the modes.

    """

    def __init__(self, femSpace, mesh):
        """
        Initializes the ModeCollection instance.

        Parameters
        ----------
        femSpace : object
            The finite element space associated with the modes.
        mesh : object
            The mesh associated with the modes.
        """
        
        self.modeList = []

        self.femSpace = femSpace
        self.mesh     = mesh
        self.names    = names


    def appendMode(self, mode):
        """
        Append an existing Mode object to the collection.

        Parameters
        ----------
        mode : Mode
            The Mode object to append.
        """

        self.modeList.append(mode)

    def appendModeFromVector(self, vector, gain=None, eigenValue=None, guess=None, waveNumber=None, frequency=None, isAdjoint=False):
        """
        Create a Mode from a coefficient vector and properties, and append it to the collection.

        Parameters
        ----------
        vector : array-like
            Coefficient vector for the mode.
        gain : float, optional
            Gain associated with the mode.
        eigenValue : complex, optional
            Eigenvalue associated with the mode.
        guess : any, optional
            Initial guess or parameter for the mode.
        waveNumber : float, optional
            Wave number associated with the mode.
        frequency : float, optional
            Frequency associated with the mode.
        isAdjoint : bool, optional
            Whether the mode is an adjoint mode (default is False).
        """

        mode = Mode(self.femSpace, self.mesh, name=self.names, isStateVector = True)
        #mode.setCoefficientArray(vector)
        mode.function.x.array[:] = vector

        mode.setGain(gain)
        mode.setEigenValue(eigenValue)
        mode.setGuess(guess)
        mode.setWaveNumber(waveNumber)
        mode.setFrequency(frequency)
        mode.isAdjoint = isAdjoint

        self.appendMode(mode)

    def appendSolutionOfEigenProblem(self, solution, guess, adjoint=False):
        """
        Append all modes from an eigenproblem solution to the collection.

        Parameters
        ----------
        solution : tuple
            Tuple containing eigenvalues, eigenvectors, and errors.
        guess : any
            Initial guess or parameter for the modes.
        adjoint : bool, optional
            Whether the modes are adjoint modes (default is False).
        """

        [eigVals, eigVecs, error] = solution
        numberOfModes = len(eigVals)

        for i in range(numberOfModes):
            mode = Mode(self.femSpace, self.mesh, name = self.names, isStateVector = True)
            mode.isAdjoint = adjoint
            mode.setError(error[i])
            mode.setGuess(guess)
            mode.setEigenValue(eigVals[i])
            mode.function.x.array[:] = eigVecs[i,:]
            self.modeList.append(mode)

    def getMaximumError(self):
        """
        Return the maximum error among all modes in the collection.

        Returns
        -------
        float
            The maximum error value.
        """

        import numpy as np
        error = []
        for mode in self.modeList:
            error.append(mode.getError())
        if len(self.modeList)==0:
            pass #TODO: throw error
        else:
            return np.amax(error)

    def getDirectEigenValueSpectrum(self):
        """
        Return the eigenvalue spectrum of direct (non-adjoint) modes.

        Returns
        -------
        list
            List of eigenvalues for direct modes.
        """

        spectrum = []
        for mode in self.modeList:
            if not mode.isAdjoint:
                spectrum.append(mode.getEigenValue())
        return spectrum


    def getNearestMode(self, gain=None, eigenValue=None, guess=None, waveNumber=None, frequency=None):
        """
        Find the nearest mode to the specified parameters.

        Parameters
        ----------
        gain : float, optional
            Gain to match.
        eigenValue : complex, optional
            Eigenvalue to match.
        guess : any, optional
            Guess to match.
        waveNumber : float, optional
            Wave number to match.
        frequency : float, optional
            Frequency to match.

        Returns
        -------
        Mode or None
            The nearest matching Mode object, or None if not found.

        Notes
        -----
        This method is not yet implemented.
        """

        #ToDo: get nearest mode to one of the above. Change handling of parameters
        pass


    def getLeadingMode(self, adjoint=False):
        """
        Return the leading mode based on the imaginary part of the eigenvalue.

        For direct modes, returns the mode with the largest imaginary part.
        For adjoint modes, returns the mode with the smallest imaginary part.

        Parameters
        ----------
        adjoint : bool, optional
            Whether to search among adjoint modes (default is False).

        Returns
        -------
        Mode or None
            The leading Mode object, or None if not found.
        """

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
        """
        Create legacy fluctuation solution objects from the modes.

        This is a wrapper for the old solution class and should be removed at the end of restructuring.

        Parameters
        ----------
        meanFlow : object
            Mean flow object.
        param : object
            Parameter object containing case information.
        FEMSpaces : object
            Finite element spaces.

        Returns
        -------
        list
            List of fluctuationSolutions objects.
        """

        fluctSolutObjList    = []
        fluctSolutObjListAdj = []

        if param.Case.AnalysisMode == "Modal":
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
        elif param.Case.AnalysisMode == "Input-Output":
            for mode in self.modeList:
                fluctSolutObjList.append(fluctuationSolutions(
                                        param,
                                        meanFlow,
                                        FEMSpaces,
                                        mode.getFrequency(),
                                        mode.function.x.array[:],
                                        mode.isAdjoint==False,
                                        0,
                                        mode.getGain()
                                         )
                )



        return fluctSolutObjList


    def getSize(self):
        """
        Return the number of modes in the collection.

        Returns
        -------
        int
            Number of modes in the collection.
        """

        return len(self.modeList)

    def popList(self):
        """
        Remove and return the last mode in the collection.

        Returns
        -------
        Mode
            The last Mode object in the collection.
        """

        return self.modeList.pop()
