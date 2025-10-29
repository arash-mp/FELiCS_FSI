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

    def appendModeFromVector(self, vector, gain=None, eigenValue=None, guess=None, waveNumber=None, frequency=None, isAdjoint=False, name = None, m=None):
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
        
        # Check that we are in Modal analysis
        if self.analysis != 'Modal':
            logger.error('appendSolutionOfEigenProblem called for non-Modal analysis in ModeCollection.')
            return

        [eigVals, eigVecs, error] = solution
        numberOfModes = len(eigVals)

        for i in range(numberOfModes):
            mode = Mode(
                self.femSpace, 
                self.mesh, 
                name = name, 
                isStateVector = True,   # NOTE: always True?
                analysis = 'Modal'
                )
            mode.isAdjoint = adjoint
            mode.setError(error[i])
            mode.setGuess(guess)
            mode.setEigenValue(eigVals[i])
            mode.function.x.array[:] = eigVecs[i,:]
            self.modeList.append(mode)
            


    def appendSolutionOfSVDProblem(self, forcingArray, omega, gains, resolventOperator, name=None):
        """
        Append all forcing modes from a SVD solution to the collection for one frequency.
        Then compute the corresponding response modes and append them as well.
        This assumes that we are doing a 'Resolvent' analysis.

        Parameters
        ----------
        forcingArray : array-like
            The forcing array.
        omega : float
            The frequency.
        gains : array-like
            The gains from the SVD solution.
        resolventOperator : FELiCS object
            The resolvent operator, instance of the ResolventOperator class.
        name : str, optional
            The name of the mode.
        """

        # Check that we are in Resolvent analysis
        if self.analysis != 'Resolvent':
            logger.error('appendSolutionOfSVDProblem called for non-Resolvent analysis in ModeCollection.')
            return
        
        # Get required operators from the resolvent operator
        if not hasattr(resolventOperator, '_W_forcing'):
            logger.error('Resolvent operator missing _W_forcing attribute in _computeResolventResponseFromForcing.')
            return
        else:
            W_forcing                       = resolventOperator._W_forcing
        if not hasattr(resolventOperator, '_W_FEM'):
            logger.error('Resolvent operator missing _W_FEM attribute in _computeResolventResponseFromForcing.')
            return
        else:
            W_FEM                           = resolventOperator._W_FEM
        if not hasattr(resolventOperator, '_P_forcing'):
            logger.error('Resolvent operator missing _P_forcing attribute in _computeResolventResponseFromForcing.')
            return
        else:
            P_forcing                       = resolventOperator._P_forcing
        if not hasattr(resolventOperator, 'getKSP'):
            logger.error('Resolvent operator missing getKSP method in _computeResolventResponseFromForcing.')
            return
        
        # Number of solutions
        numberOfSolutions                   = gains.shape[0]
        
        # For each forcing, append to mode and compute response
        logger.debug(f'Appending {numberOfSolutions} forcing and response modes for omega={omega}.')
        for i in range(numberOfSolutions):
            # Get petsc vectors from petsc matrices (=get petsc vectors with correct sizes)             
            X1, X2                          = W_forcing.getVecs()             
            X1.setValues(range(0,len(forcingArray)), forcingArray[:,i])             
            Y1, Y2                          = W_FEM.getVecs()             

            # Solve forcings = Pu*eigenVectors             
            P_forcing.mult(X1,Y1)             
            forcings                        = Y1.getValues(range(0, Y1.getSize()))

            # Setting forcing into mode object
            modeForcing                     = Mode(
                self.femSpace, 
                self.mesh, 
                name                        = name, 
                isStateVector               = True, # NOTE: always True?
                analysis                    = 'Resolvent'
                )
            modeForcing.setGain(np.real(gains[i]))  # NOTE: These are gains squared
            modeForcing.setGainNumber(i)
            modeForcing.setFrequency(omega)
            modeForcing.function.x.array[:] = forcings
            self.modeList.append(modeForcing)
            
            # Compute response
            # Solve Y1 = -1j * B_femWeight * forcings
            W_FEM.mult(Y1,Y2)             
            Y2.scale(-1j)             
            # Solve (A-omega*B)*responses = Y1
            resolventOperator.getKSP().solve(Y2,Y1)
            responses                       = Y1.getValues(range(0, Y1.getSize()))
            
            # Setting response into mode object
            modeForcing                     = Mode(
                self.femSpace, 
                self.mesh, 
                name                        = name, 
                isStateVector               = True,
                analysis                    = 'Resolvent'
                )
            modeForcing.setGain(np.real(gains[i]))  # NOTE: These are gains squared
            modeForcing.setGainNumber(i)
            modeForcing.setFrequency(omega)
            modeForcing.function.x.array[:] = responses
            modeForcing.isResponse          = True
            self.modeList.append(modeForcing)


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
    
    def getSpectrum(self):
        """
        Returns the spectrum for either `Modal` or `Resolvent` type of analysis.
        
        
        Returns:
        -------
        spectrum: array-like 
            Array of eigenvalues or gains at each omega and either (i) Direct/Adjoint for 'Modal'
            or (ii) GainNumber for `Resolvent`.
        header: list
            List of headers corresponding to the spectrum values.
        """
        
        # Mode list
        modeList        = self.modeList
        
        # For modal analysis
        if self.analysis == 'Modal':
            hasAdjoint  = any(mode.isAdjoint for mode in modeList)
            nLines      = len(modeList)//2 if hasAdjoint else len(modeList)
            nCols       = 4 if hasAdjoint else 2
            
            # Define the header
            if hasAdjoint:
                header  = ['omega_direct_r','omega_direct_i','omega_adjoint_r','omega_adjoint_i']
            else:
                header  = ['omega_direct_r','omega_direct_i']
                
            # Define the spectrum array
            spectrum    = np.zeros((nLines, nCols), dtype=float)
            ctr_line_direct     = 0
            ctr_line_adjoint    = 0
            for i, mode in enumerate(modeList):
                if not mode.isAdjoint:
                    spectrum[ctr_line_direct,0] = mode.getEigenValue().real
                    spectrum[ctr_line_direct,1] = mode.getEigenValue().imag
                    ctr_line_direct += 1
                else:
                    spectrum[ctr_line_adjoint,2] = mode.getEigenValue().real
                    spectrum[ctr_line_adjoint,3] = mode.getEigenValue().imag
                    ctr_line_adjoint += 1
            
        # Resolvent or IO case (consider only the response modes)
        elif self.analysis in ['Resolvent', 'Input-Output']:
            Ncols           = 2 + max(mode.getGainNumber() for mode in modeList) # Mode numbers start at 0
            frequencyList   = np.unique([mode.getFrequency() for mode in modeList])
            NLines          = len(frequencyList)

            # Define the header
            header          = ["omega"]
            header.extend([f'gain_{i}' for i in range(Ncols - 1)])

            # Define the spectrum array
            spectrum    = np.zeros((NLines, Ncols), dtype=float)
            hasResponse = any(mode.isResponse for mode in modeList)
            # Only loop on one type of modes to avoid duplicates
            modeTypeToLoop = 'Response' if hasResponse else 'Forcing'
            for i, mode in enumerate(modeList):
                if (modeTypeToLoop == 'Response' and mode.isResponse) or (modeTypeToLoop == 'Forcing' and not mode.isResponse):
                    currentFrequency    = mode.getFrequency()
                    LineIndex           = np.where(frequencyList == currentFrequency)[0]
                    spectrum[LineIndex[0], 0] = currentFrequency
                    spectrum[LineIndex[0], 1 + mode.getGainNumber()] = mode.getGain()

        return spectrum, header

    def exportSpectrumToCSV(self, filePath):
        """
        Export the spectrum to a CSV file.

        Parameters
        ----------
        filePath : str
            Path to the output CSV file.
        """

        spectrum, header = self.getSpectrum()
        np.savetxt(
            filePath, 
            spectrum, 
            delimiter=',', 
            header=','.join(header), 
            comments=''
        )
        if self.analysis == 'Modal':
            logger.info(f'Eigenvalue spectrum exported to {filePath}.')
        else:
            logger.info(f'Gains spectrum exported to {filePath}.')

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
        elif param.Case.AnalysisMode == "Resolvent":
            for mode in self.modeList:
                fluctSolutObjList.append(fluctuationSolutions(                                 
                    param,                                 
                    meanFlow,                                 
                    FEMSpaces,                                 
                    mode.getFrequency(),                                 
                    mode.function.x.array[:],                                
                    False if not mode.isResponse else True,                                 
                    mode.getGainNumber(),                                 
                    mode.getGain(),                                 
                ))
            
        else:
            logger.error("Analysis type not recognized in getOldSolutionObject of ModeCollection.")
            raise RuntimeError("Analysis type not recognized in getOldSolutionObject of ModeCollection.")



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

        TODO: return the gainNumber for resolvent case (for IO it's always 0)
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
        gainNumbersModeFiles    = []
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
                # Get the gain number from filename
                gainNumberStr       = f.split('gain')[1].split('.h5')[0]
                gainNumbersModeFiles.append(int(gainNumberStr))

            elif self.analysis == 'Input-Output':
                gainNumbersModeFiles.append(0)  # Always 0 for IO modes
                if responsePattern in f:
                    typesModeFiles.append('Response')
                else:
                    logger.error(f'Could not determine mode type from filename "{f}". Skipping this file.')
                    raise RuntimeError(f'Could not determine mode type from filename "{f}".')
            
            # Extract omega from filename
            # NOTE: different structure of name for Modal and Resolvent TODO: unify?
            if self.analysis == 'Modal':
                omegaStr            = f.split(filePrefix+typesModeFiles[i]+'_')[1].split('.h5')[0]
            else:
                omegaStr            = f.split(filePrefix)[1].split('_')[0]
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

        return modeFiles, omegasModeFiles, typesModeFiles, gainNumbersModeFiles

    def importData(self, reader, importFolder, omegas=None, modeType=None, gainNumber=None):
        """
        Import mode collection data from a specified folder using a reader.

        Parameters
        ----------
        reader : Reader
            The Reader object used to read the mode data.
        importFolder : str
            Path to the folder containing the mode files.
        omegas : list of float, optional
            List of frequencies to import. If None, all frequencies are imported.
        modeType : str, optional
            Type of mode to import (e.g., 'Direct', 'Adjoint', 'Response', 'Forcing'). If None, all types are imported.
        gainNumber : int, optional
            Gain number to import (for Resolvent analysis). If None, all gain numbers are imported
        """
        
        # Get the list of mode files in the directory, omegas values, and mode types
        h5Files, fileOmegas, fileTypes, gainNb  = self._getAndSortModeFilesInDir(importFolder)
        
        # If omega was given as input, filter files accordingly
        if omegas is not None:
            matchingOmegasIndices       = [i for i, omega in enumerate(fileOmegas) if omega in np.round(omegas, 3)]
            if len(matchingOmegasIndices) == 0:
                logger.error('No matching omegas found in the import folder for the specified omegas.')
                return
            h5Files                     = [h5Files[i] for i in matchingOmegasIndices]
            logger.info(f'Importing only modes for omegas: {np.round([fileOmegas[i] for i in matchingOmegasIndices],3)}')

        # If a mode type was given as input, filter files accordingly
        if modeType is not None:
            matchingTypeIndices         = [i for i, mType in enumerate(fileTypes) if mType == modeType]
            if len(matchingTypeIndices) == 0:
                logger.error(f'No matching mode types found in the import folder for the specified type "{modeType}".')
                return
            h5Files                     = [h5Files[i] for i in matchingTypeIndices]
            logger.info(f'Importing only modes with type: {modeType}')
            
        # If a gain number was given as input, filter files accordingly (for Resolvent analysis)
        if gainNumber is not None and self.analysis == 'Resolvent':
            matchingGainIndices         = [i for i, gNum in enumerate(gainNb) if gNum == gainNumber]
            if len(matchingGainIndices) == 0:
                logger.error(f'No matching gain numbers found in the import folder for the specified gain number "{gainNumber}".')
                return
            h5Files                     = [h5Files[i] for i in matchingGainIndices]
            logger.info(f'Importing only modes with gain number: {gainNumber}')

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
            
            # For resolvent an IO, set gainNumber before importing
            if self.analysis in ['Resolvent', 'Input-Output']:
                mode.setGainNumber(gainNb[i])
                
            # For resolvent, set isResponse based on filename
            # TODO: use a more general "modeType" attribute in Mode?
            if self.analysis == 'Resolvent':
                if 'Response' in h5Files[i]:
                    mode.isResponse = True
                else:
                    mode.isResponse = False
            elif self.analysis == 'Input-Output':
                mode.isResponse = True  # Always response for IO modes
            
            mode.importData(
                reader,
                importFolder,
                importFile = h5Files[i],
            )
            logger.info(f'Imported mode {i+1}/{numberOfModes} from "{h5Files[i]}".')
            self.modeList.append(mode)