import  os
import  numpy                   as np
from    .Mode                   import Mode
from    .fluctuationClass       import fluctuationSolutions
from 	FELiCS.Misc.logging     import Logger

# Get the logger
logger = Logger.get_logger("felics")

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

    def __init__(self, femSpace, mesh, isStateVector=True, m=0, analysis='Modal'):
        """
        Initializes the ModeCollection instance.

        Parameters
        ----------
        femSpace : object
            The finite element space associated with the modes.
        mesh : object
            The mesh associated with the modes.
        TODO: complete docstring
        """

        self.modeList       = []
        self.femSpace       = femSpace
        self.mesh           = mesh
        self.isStateVector  = isStateVector
        self.m              = m
        self.analysis       = analysis


    def appendMode(self, mode):
        """
        Append an existing Mode object to the collection.

        Parameters
        ----------
        mode : Mode
            The Mode object to append.
        """

        self.modeList.append(mode)

    def appendModeFromVector(self, vector, gain=None, eigenValue=None, guess=None, waveNumber=None, frequency=None, isAdjoint=False, name = [], m=None):
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

        mode = Mode(self.femSpace, self.mesh, name, isStateVector = True)
        #mode.setCoefficientArray(vector)
        mode.function.x.array[:] = vector

        mode.setGain(gain)
        mode.setEigenValue(eigenValue)
        mode.setGuess(guess)
        mode.setWaveNumber(waveNumber)
        mode.setFrequency(frequency)
        mode.isAdjoint = isAdjoint

        self.appendMode(mode)

    def appendSolutionOfEigenProblem(self, solution, guess, adjoint=False, name=None):
        """
        Append all modes from an eigenproblem solution to the collection.
        This assumes that we are doing a 'Modal' analysis.

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
            mode = Mode(
                self.femSpace, 
                self.mesh, 
                name = name, 
                isStateVector = True,
                analysis = 'Modal'
                )
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

    def _getAndSortModeFilesInDir(self, importFolder):
        """
        Get and sort mode files in a specified directory.

        Parameters
        ----------
        importFolder : str
            Path to the folder containing the mode files.

        Returns
        -------
        list
            Sorted list of mode file names.
        """

        # List of all h5 files in the folder
        h5Files             = [f for f in os.listdir(importFolder) if f.endswith('.h5')]
        
        # File patterns based on analysis type
        if self.analysis == 'Modal':
            filePrefix      = 'ModalSolution_Omega_'
            directPattern   = 'Direct'
            adjointPattern  = 'Adjoint'
        elif self.analysis == 'Resolvent':
            filePrefix      = 'Resolvent_Omega'
            responsePattern = 'Response'
            forcingPattern  = 'Forcing'
        elif self.analysis == 'Input-Output':
            filePrefix      = 'Input-Output_Omega'
            responsePattern = 'Response'
        else:
            logger.error(f'Analysis type "{self.analysis}" not recognized. Cannot import mode collection.')
            return
        
        # List of files that match the prefix
        modeFiles               = [f for f in h5Files if filePrefix in f]
        
        # Go over files and classify them
        omegasModeFiles         = []
        typesModeFiles          = []
        for i, f in enumerate(modeFiles):
            
            # First we get the type
            if self.analysis == 'Modal':
                if directPattern in f:
                    typesModeFiles.append('Direct')
                elif adjointPattern in f:
                    typesModeFiles.append('Adjoint')
                else:
                    logger.error(f'Could not determine mode type from filename "{f}".')
                    raise RuntimeError(f'Could not determine mode type from filename "{f}".')
            elif self.analysis == 'Resolvent':
                if responsePattern in f:
                    typesModeFiles.append('Response')
                elif forcingPattern in f:
                    typesModeFiles.append('Forcing')
                else:
                    logger.error(f'Could not determine mode type from filename "{f}". Skipping this file.')
                    raise RuntimeError(f'Could not determine mode type from filename "{f}".')
            elif self.analysis == 'Input-Output':
                if responsePattern in f:
                    typesModeFiles.append('Response')
                else:
                    logger.error(f'Could not determine mode type from filename "{f}". Skipping this file.')
                    raise RuntimeError(f'Could not determine mode type from filename "{f}".')
            
            # Extract omega from filename
            omegaStr            = f.split(filePrefix+typesModeFiles[i]+'_')[1].split('.h5')[0]
            # Can be float (no "j") or a complex number (with "j")
            if 'j' in omegaStr:
                try:
                    omega       = complex(omegaStr)
                except Exception as e:
                    logger.error(f'Could not parse omega from filename "{f}". Skipping this file. Error: {e}')
                    continue
            else:
                try:
                    omega       = float(omegaStr)
                except Exception as e:
                    logger.error(f'Could not parse omega from filename "{f}". Skipping this file. Error: {e}')
                    continue
            omegasModeFiles.append(omega)
        
        return modeFiles, omegasModeFiles, typesModeFiles

    def importData(self, reader, importFolder, omegas=None, modeType=None):
        """
        Import mode collection data from a specified folder using a reader.

        Parameters
        ----------
        reader : object
            Reader object to handle data import.
        importFolder : str
            Path to the folder containing the mode collection data.
        omegas : list, optional
            List of eigenfrequencies to import. If None, all modes are imported.
        """
        
        # Get the list of mode files in the directory, omegas values, and mode types
        h5Files, fileOmegas, fileTypes  = self._getAndSortModeFilesInDir(importFolder)
        
        # If omega was given as input, filter files accordingly
        if omegas is not None:
            matchingOmegasIndices       = [i for i, omega in enumerate(fileOmegas) if omega in np.round(omegas, 3)]
            if len(matchingOmegasIndices) == 0:
                logger.error('No matching omegas found in the import folder for the specified omegas.')
                return
            h5Files                     = [h5Files[i] for i in matchingOmegasIndices]

        # If a mode type was given as input, filter files accordingly
        if modeType is not None:
            matchingTypeIndices         = [i for i, mType in enumerate(fileTypes) if mType == modeType]
            if len(matchingTypeIndices) == 0:
                logger.error(f'No matching mode types found in the import folder for the specified type "{modeType}".')
                return
            h5Files                     = [h5Files[i] for i in matchingTypeIndices]

        # Import each mode file
        numberOfModes                   = len(h5Files)
        logger.info(f'Importing {numberOfModes} modes into mode collection.')
        for i in range(numberOfModes):
            mode                        = Mode(
                self.femSpace, 
                self.mesh, 
                isStateVector=self.isStateVector, 
                m=self.m, 
                analysis=self.analysis
            )
            mode.importData(
                reader,
                importFolder,
                importFile = h5Files[i],
            )
            logger.info(f'Imported mode {i+1}/{numberOfModes} from "{h5Files[i]}".')
            self.modeList.append(mode)