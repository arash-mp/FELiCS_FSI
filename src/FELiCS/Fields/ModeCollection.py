#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
# Standard libraries
import  os

# Third party libraries
import  h5py

import numpy as np

# Local Libraries and methods
from    FELiCS.Fields.FluctuationClass  import FluctuationSolutions
from    FELiCS.Fields.Mode              import (
    Mode,
    AnalysisType,
    ModeType,
)
from 	FELiCS.Misc.logging             import Logger,  log_and_raise


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
    femSpace : dolfinx.fem.FunctionSpace
        The finite element space associated with the modes.
    mesh : FELiCS.SpaceDisc.FELiCSMesh
        The mesh associated with the modes.
    isStateVector : bool, optional
        If True, modes in this collection are treated as state vectors in a
        mixed space. Default is True.
    analysisType : str, optional
        Type of analysis this collection belongs to. The string is mapped
        (case-insensitively) to :class:`AnalysisType` (e.g. 'Modal', 'Resolvent',
        'Input_Output'). Default is ``'Modal'``.

    Attributes
    ----------
    modeList : list
        List of Mode objects in the collection.
    femSpace : dolfinx.fem.FunctionSpace
        The finite element space associated with the modes.
    mesh : FELiCS.SpaceDisc.FELiCSMesh
        The mesh associated with the modes.
    isStateVector : bool
        Indicates if modes are treated as state vectors.
    analysisType : AnalysisType
        The specific analysis type enum associated with this collection.
    """

    def __init__(
        self,
        femSpace,
        mesh,
        isStateVector=True,
        analysisType='Modal'
    ):
        """
        Initialize a ModeCollection instance.

        Parameters
        ----------
        femSpace : dolfinx.fem.FunctionSpace
            Finite element space associated with the modes.
        mesh : FELiCS.SpaceDisc.FELiCSMesh
            Mesh on which the modes are defined.
        isStateVector : bool, optional
            If True, modes in this collection are treated as state vectors in a
            mixed space. Default is True.
        analysisType : {'Modal', 'Resolvent', 'Input_Output'}, optional
            Type of analysis this collection belongs to. The string is mapped
            (case-insensitively) to :class:`AnalysisType`. Default is ``'Modal'``.
        """

        self.modeList       = []
        self.femSpace       = femSpace
        self.mesh           = mesh
        self.isStateVector  = isStateVector
        self.analysisType   = AnalysisType[analysisType.upper()]


    def describe(
        self
    ):
        """
        Print a description of the ModeCollection.
        """

        logger.info(f'ModeCollection for {self.analysisType} analysis.')
        logger.info(f'Contains {len(self.modeList)} modes. Contents:')
        
        # Loop over modes and print their description
        for i, mode in enumerate(self.modeList):
            logger.info(f'  Mode {i}:')
            mode.describe()

    def append_mode(
        self,
        mode
    ):
        """
        Append an existing Mode object to the collection.

        Parameters
        ----------
        mode : Mode
            The Mode object to append.
        """

        self.modeList.append(mode)

    def append_mode_from_vector(
        self,
        vector,
        gain=None,
        eigen_value=None,
        guess=None,
        wave_number=None,
        frequency=None,
        modeType = None,
        m=0
    ):
        """
        Create a Mode from a coefficient vector and properties, and append it
        to the collection.

        Parameters
        ----------
        vector : array-like
            Coefficient vector for the mode.
        gain : float, optional
            Gain associated with the mode (mainly for resolvent analysis).
        eigenValue : complex, optional
            Eigenvalue associated with the mode (for modal analysis).
        guess : any, optional
            Initial guess or parameter used to obtain this mode.
        waveNumber : float, optional
            Wave number associated with the mode.
        frequency : float, optional
            Frequency associated with the mode (for resolvent / input–output).
        modeType : str or None, optional
            Mode role within the chosen analysis type, mapped to
            :class:`ModeType` (e.g. ``'direct'``, ``'adjoint'``,
            ``'response'``, ``'forcing'``). If ``None``, a default is chosen
            inside :class:`Mode` based on the analysis type.
        m : int, optional
            Azimuthal wavenumber associated with the mode. Default is 0.
        """
        mode = Mode(
            self.femSpace,
            self.mesh,
            isStateVector = True,
            analysisType = self.analysisType.name,
            modeType = modeType,
            m = m,
        )
        mode.set_coefficient_array(vector)

        mode.gain           = gain
        mode.eigen_value    = eigen_value
        mode.guess          = guess
        mode.wave_number    = wave_number
        mode.frequency      = frequency

        self.append_mode(mode)

    def append_solution_of_eigen_problem(
        self,
        solution,
        guess,
        adjoint=False,
        name=None,
        m=0
    ):
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
        name : str, optional
            A custom name for the mode.
        m : int, optional
            Azimuthal wavenumber associated with the mode. Default is 0.
        """

        # Check that we are in Modal analysis
        if self.analysisType != AnalysisType.MODAL:
            logger.error('appendSolutionOfEigenProblem called for non-Modal analysis in ModeCollection.')
            return

        [eigVals, eigVecs, error] = solution
        numberOfModes   = len(eigVals)

        if not adjoint:
            modeType    = 'direct'
        else:
            modeType    = 'adjoint'
            guess       = np.conj(guess)

        for i in range(numberOfModes):
            mode = Mode(
                self.femSpace,
                self.mesh,
                name = name,
                isStateVector = True,
                analysisType = 'Modal',
                modeType = modeType,
                m = m
            )
            mode.isAdjoint      = adjoint
            mode.error          = error[i]
            mode.guess          = guess
            mode.eigen_value    = eigVals[i]
            mode.set_coefficient_array(eigVecs[i,:])
            self.modeList.append(mode)

    def append_solution_of_svd_problem(
        self,
        forcingArray,
        omega,
        gains,
        resolventOperator,
        name=None,
        m=0
    ):
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
        m : int, optional
            Azimuthal wavenumber associated with the mode. Default is 0.
        """

        # Check that we are in Resolvent analysis
        if self.analysisType != AnalysisType.RESOLVENT:
            logger.error('appendSolutionOfSVDProblem called for non-Resolvent analysis in ModeCollection.')
            return
        
        # Fetch required operators from resolventOperator
        W_forcing                       = getattr(
            resolventOperator,
            '_W_forcing',
            None,
        )
        W_FEM                           = getattr(
            resolventOperator,
            '_W_FEM',
            None,
        )
        P_forcing                       = getattr(
            resolventOperator,
            '_P_forcing',
            None,
        )
        getKSP                          = getattr(
            resolventOperator,
            'getKSP',
            None,
        )

        if W_forcing is None or W_FEM is None or P_forcing is None or getKSP is None:
            missing = []
            if W_forcing is None: 
                missing.append('_W_forcing')
            if W_FEM is None:     
                missing.append('_W_FEM')
            if P_forcing is None: 
                missing.append('_P_forcing')
            if getKSP is None:    
                missing.append('getKSP')
            logger.error(f'Resolvent operator missing required attributes/methods: {", ".join(missing)}.')
            return

        # Number of solutions
        numberOfSolutions               = gains.shape[0]
        
        # For each forcing, append to mode and compute response
        logger.debug(f'Appending {numberOfSolutions} forcing and response modes for omega={omega}.')
        for i in range(numberOfSolutions):

            # Get petsc vectors from petsc matrices (=get petsc vectors with correct sizes)
            X1, X2                      = W_forcing.getVecs()             
            X1.setValues(
                range(
                    0,
                    len(forcingArray)
                ),
                forcingArray[:,i]
            )             
            Y1, Y2                      = W_FEM.getVecs()             

            # Solve forcings = Pu*eigenVectors             
            P_forcing.mult(
                X1,
                Y1
            )             
            forcings                    = Y1.getValues(range(
                0,
                Y1.getSize(),
            ))

            # Setting forcing into mode object
            modeForcing                 = Mode(
                self.femSpace,
                self.mesh,
                name                    = name,
                isStateVector           = True,
                analysisType            = 'Resolvent',
                modeType                = 'forcing',
                m                       = m,
            )
            modeForcing.gain            = np.real(gains[i])  # NOTE: These are gains squared
            modeForcing.gain_number     = i
            modeForcing.frequency       = omega
            modeForcing.set_coefficient_array(forcings)
            self.modeList.append(modeForcing)
            
            # Compute response
            # Solve Y1 = -1j * B_femWeight * forcings
            W_FEM.mult(
                Y1,
                Y2
            )             
            Y2.scale(-1j)             
            # Solve (A-omega*B)*responses = Y1
            resolventOperator.getKSP().solve(
                Y2,
                Y1
            )
            responses                   = Y1.getValues(range(
                0,
                Y1.getSize()
            ))
            
            # Setting response into mode object
            modeResponse                = Mode(
                self.femSpace,
                self.mesh,
                name                    = name,
                isStateVector           = True,
                analysisType            = 'Resolvent',
                modeType                = 'response',
                m                       = m
            )
            modeResponse.gain           = np.real(gains[i])  # NOTE: These are gains squared
            modeResponse.gain_number    = i
            modeResponse.frequency      = omega
            modeResponse.set_coefficient_array(responses)
            self.modeList.append(modeResponse)


    def get_maximum_error(
        self
    ):
        """
        Return the maximum error among all modes in the collection.

        Returns
        -------
        float
            The maximum error value.
        """
        error = []
        for mode in self.modeList:
            error.append(mode.error)
        if len(self.modeList)==0:
            pass #TODO: throw error
        else:
            return np.amax(error)

    def get_direct_eigen_value_spectrum(
        self
    ):
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
                spectrum.append(mode.eigen_value)
        return spectrum
    
    def get_spectrum(
        self
    ):
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

        spectrum = []
        header   = []

        # Mode list
        modeList = self.modeList

        # For modal analysis
        if self.analysisType == AnalysisType.MODAL:
            hasAdjoint  = any(mode.modeType == ModeType.ADJOINT for mode in modeList)
            nLines      = len(modeList)//2 if hasAdjoint else len(modeList)
            nCols       = 4 if hasAdjoint else 2
            
            # Define the header
            if hasAdjoint:
                header  = ['omega_direct_r','omega_direct_i','omega_adjoint_r','omega_adjoint_i']
            else:
                header  = ['omega_direct_r','omega_direct_i']
                
            # Define the spectrum array
            spectrum    = np.zeros(
                (nLines, nCols),
                dtype=float
            )
            ctr_line_direct     = 0
            ctr_line_adjoint    = 0
            for i, mode in enumerate(modeList):
                if not mode.modeType == ModeType.ADJOINT:
                    spectrum[ctr_line_direct,0] = mode.eigen_value.real
                    spectrum[ctr_line_direct,1] = mode.eigen_value.imag
                    ctr_line_direct += 1
                else:
                    spectrum[ctr_line_adjoint,2] = mode.eigen_value.real
                    spectrum[ctr_line_adjoint,3] = mode.eigen_value.imag
                    ctr_line_adjoint += 1
            
        # Resolvent or IO case (consider only the response modes)
        elif self.analysisType in [AnalysisType.RESOLVENT]:
            Ncols           = 2 + max(mode.gain_number for mode in modeList) # Mode numbers start at 0
            frequencyList   = np.unique([mode.frequency for mode in modeList])
            NLines          = len(frequencyList)

            # Define the header
            header          = ["omega"]
            header.extend([f'gain_{i}' for i in range(Ncols - 1)])

            # Define the spectrum array
            spectrum    = np.zeros(
                (NLines, Ncols),
                dtype=float,
            )
            hasResponse = any(mode.modeType == ModeType.RESPONSE for mode in modeList)
            
            # Only loop on one type of modes to avoid duplicates
            modeTypeToLoop = 'Response' if hasResponse else 'Forcing'
            for i, mode in enumerate(modeList):
                if (modeTypeToLoop == 'Response' and mode.modeType == ModeType.RESPONSE) or (modeTypeToLoop == 'Forcing' and mode.modeType == ModeType.FORCING):
                    currentFrequency    = mode.frequency
                    LineIndex           = np.where(frequencyList == currentFrequency)[0]
                    spectrum[LineIndex[0], 0] = currentFrequency
                    spectrum[LineIndex[0], 1 + mode.gain_number] = mode.gain

        return spectrum, header

    def export_spectrum_to_csv(
        self,
        writer
    ):
        """
        Export the spectrum to a CSV file.

        Parameters
        ----------
        writer: FELiCS.IO.Writer object
        """

        if self.analysisType == AnalysisType.MODAL:
            filePath = os.path.join(
                writer.exportFolder,
                "spectrum.csv"
            )
        else:
            filePath = os.path.join(
                writer.exportFolder,
                "gains.csv"
            )

        spectrum, header = self.get_spectrum()
        np.savetxt(
            filePath,
            spectrum,
            delimiter=',',
            header=','.join(header),
            comments=''
        )
        if self.analysisType == AnalysisType.MODAL:
            logger.info(f'Eigenvalue spectrum exported to {filePath}.')
        else:
            logger.info(f'Gains spectrum exported to {filePath}.')

    def get_nearest_mode(
        self,
        gain=None,
        eigen_value=None,
        guess=None,
        wave_number=None,
        frequency=None
    ):
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

        Notes
        -----
        This method is not yet implemented and will be added in future releases.
        """

        #TODO: get nearest mode to one of the above. Change handling of parameters
        #TODO: until then: throw ERROR!!!!
        pass


    def get_leading_mode(
        self,
        adjoint=False
    ):
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
                    eigen_value = mode.eigen_value
                    if np.imag(eigen_value) > growthRateMax:
                        growthRateMax = np.imag(eigen_value)
                        leadingMode = mode
        elif adjoint:
            growthRateMin = 9990.
            for mode in self.modeList:
                if mode.isAdjoint == adjoint:
                    eigen_value = mode.eigen_value
                    if np.imag(eigen_value) < growthRateMin:
                        growthRateMin = np.imag(eigen_value)
                        leadingMode = mode
                    
        return leadingMode


    def get_old_solution_object(
        self,
        meanFlow,
        param,
        FEMSpaces
    ):
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

        Notes
        -------
        In future relases, this method will be deleted
        """

        fluctSolutObjList    = []
        fluctSolutObjListAdj = []

        if param.Case.AnalysisMode == "Modal":
            for mode in self.modeList:
                # construct for each EVal and EVec a fluctuationSolution
                fluctSolutObjList.append(FluctuationSolutions(
                    param,
                    meanFlow,
                    FEMSpaces,
                    mode.eigen_value,
                    mode.function.x.array[:],
                    mode.isAdjoint==False,
                ))
        elif param.Case.AnalysisMode == "Input-Output":
            for mode in self.modeList:
                fluctSolutObjList.append(FluctuationSolutions(
                    param,
                    meanFlow,
                    FEMSpaces,
                    mode.frequency,
                    mode.function.x.array[:],
                    mode.isAdjoint==False,
                    0,
                    mode.gain,
                ))
        elif param.Case.AnalysisMode == "Resolvent":
            for mode in self.modeList:
                fluctSolutObjList.append(FluctuationSolutions(
                    param,
                    meanFlow,
                    FEMSpaces,
                    mode.frequency,
                    mode.function.x.array[:],
                    False if not mode.isResponse else True,
                    mode.gain_number,
                    mode.gain,
                ))
            
        else:
            log_and_raise(logger, "Analysis type not recognized in getOldSolutionObject of ModeCollection.", RuntimeError)
        return fluctSolutObjList


    def getSize(
        self
    ):
        """
        Return the number of modes in the collection.

        Returns
        -------
        int
            Number of modes in the collection.
        """

        return len(self.modeList)

    def pop_list(
        self
    ):
        """
        Remove and return the last mode in the collection.

        Returns
        -------
        Mode
            The last Mode object in the collection.
        """

        return self.modeList.pop()

    def _get_and_sort_mode_files_in_dir(
        self,
        importFolder
    ):
        """
        Get and classify mode files in a specified directory.

        Parameters
        ----------
        importFolder : str
            Path to the folder containing the mode files.

        Returns
        -------
        modeFiles : list of str
            List of HDF5 file names that match the expected naming pattern for
            the current analysis type.
        omegasModes : list of float or complex
            List of omega values parsed from the corresponding file names.
        typesModeFiles : list of str
            List of mode type strings (e.g. 'Direct', 'Adjoint',
            'Response', 'Forcing') inferred from the file names.
        gainNumbersModeFiles : list of int
            List of gain numbers parsed from the file names (0 for modal and
            input–output cases, or the gain index for resolvent modes).

        Notes
        -----
        For resolvent and input–output cases, the filename convention is used
        to extract both omega and gain numbers. Errors in parsing result in
        logged messages and skipping the affected files.
        """
        # List of all h5 files in the folder
        h5Files             = [f for f in os.listdir(importFolder) if f.endswith('.h5')]
        
        # File patterns based on analysis type
        if self.analysisType == AnalysisType.MODAL:
            filePrefix      = 'Mode_Modal_'
            directPattern   = 'Direct'
            adjointPattern  = 'Adjoint'
        elif self.analysisType == AnalysisType.RESOLVENT:
            filePrefix      = 'Mode_Resolvent_'
            responsePattern = 'Response'
            forcingPattern  = 'Forcing'
        elif self.analysisType == AnalysisType.INPUT_OUTPUT:
            filePrefix      = 'Mode_Input_output_'
            responsePattern = 'Response'
        else:
            logger.error(f'Analysis type "{self.analysisType}" not recognized. Cannot import mode collection.')
            return
        
        # List of files that match the prefix
        modeFiles               = [f for f in h5Files if filePrefix in f]
        
        # Go over files and classify them
        omegasModes             = []
        typesModeFiles          = []
        gainNumbersModeFiles    = []
        
        for i, f in enumerate(modeFiles): 
            if self.analysisType == AnalysisType.MODAL:
                if directPattern in f:
                    typesModeFiles.append('Direct')
                elif adjointPattern in f:
                    typesModeFiles.append('Adjoint')
                else:
                    log_and_raise(logger, f'Could not determine mode type from filename "{f}".', RuntimeError)
                # Get the complex eigenvalue from file
                with h5py.File(
                    os.path.join(
                        importFolder,
                        f
                    ),
                    'r',
                ) as file:
                    omegasModes.append(file["omega"][()])
                
            elif self.analysisType == AnalysisType.RESOLVENT:
                if responsePattern in f:
                    typesModeFiles.append('Response')
                elif forcingPattern in f:
                    typesModeFiles.append('Forcing')
                else:
                    log_and_raise(logger, f'Could not determine mode type from filename "{f}". Skipping this file.', RuntimeError)
                # Get the gain number and frequency from file
                with h5py.File(
                    os.path.join(
                        importFolder,
                        f
                    ),
                    'r',
                ) as file:
                    omegasModes.append(file["omega"][()])
                    gainNumbersModeFiles.append(int(file["number"][()]))

            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                gainNumbersModeFiles.append(0)  # Always 0 for IO modes
                if responsePattern in f:
                    typesModeFiles.append('Response')
                else:
                    log_and_raise(logger, f'Could not determine mode type from filename "{f}". Skipping this file.', RuntimeError)
                # Get the gain number and frequency from file
                with h5py.File(
                    os.path.join(
                        importFolder,
                        f,
                    ),
                    'r',
                ) as file:
                    omegasModes.append(file["omega"][()])
                    gainNumbersModeFiles.append(0) # Always 0 for IO modes

        return modeFiles, omegasModes, typesModeFiles, gainNumbersModeFiles


    def export_modes(
        self,
        writer,
        onlyNewN = 0
    ):
        """
        Export modes in the collection to HDF5 files.

        Parameters
        ----------
        writer : FELiCS.IO.Writer
            Writer object providing ``exportFieldToH5``.
        onlyNewN : int, optional
            If 0 (default), export all modes in the collection. If positive,
            export only the last ``onlyNewN`` modes that were added.
        """
        if onlyNewN == 0:  # export all modes
            start = 0  
        else:              # export only the newest modes, number given by 'onlyNewN'
            start = len(self.modeList) - onlyNewN 

        for i in range(
            start,
            len(self.modeList),
        ):
            logger.debug(f"Exporting mode {i+1}/{len(self.modeList)} to file.")
            self.modeList[i].export_to_h5(writer)


    def import_data(
        self,
        reader,
        importFolder,
        omegas=None,
        modeType=None,
        gain_number=None
    ):
        """
        Import mode collection data from a specified folder using a reader.

        Parameters
        ----------
        reader : Reader
            The Reader object used to read the mode data.
        importFolder : str
            Path to the folder containing the mode files.
        omegas : list of float, optional
            List of frequencies to import. If None, all frequencies found in
            the folder are considered (matching is done on rounded values).
        modeType : str, optional
            Type of mode to import (e.g. 'Direct', 'Adjoint', 'Response',
            'Forcing'). If None, all detected types are imported.
        gainNumber : int, optional
            Gain number to import (for Resolvent analysis). If None, all gain
            numbers are imported.
        """
        
        # Get the list of mode files in the directory, omegas values, and mode types
        h5Files, fileOmegas, fileTypes, gainNb  = self._get_and_sort_mode_files_in_dir(importFolder)
        
        # If omega was given as input, filter files accordingly
        if omegas is not None:

            #matchingOmegasIndices       = [i for i, omega in enumerate(np.round(fileOmegas, 3)) if omega in np.round(omegas, 3)]
            closestOmegasIndices        = [np.argmin(np.abs(np.array(fileOmegas) - omega)) for omega in omegas]
            if len(closestOmegasIndices) == 0:
                logger.error('No matching omegas found in the import folder for the specified omegas.')
                return
            h5Files                     = [h5Files[i] for i in closestOmegasIndices]
            logger.info(
                f'Importing modes for closest omegas: {np.round([fileOmegas[i] for i in closestOmegasIndices],3,)}'
            )
        
        # If a mode type was given as input, filter files accordingly
        if modeType is not None:
            matchingTypeIndices         = [i for i, mType in enumerate(fileTypes) if mType == modeType]
            if len(matchingTypeIndices) == 0:
                logger.error(f'No matching mode types found in the import folder for the specified type "{modeType}".')
                return
            h5Files                     = [h5Files[i] for i in matchingTypeIndices]
            logger.info(f'Importing only modes with type: {modeType}')
            
        # If a gain number was given as input, filter files accordingly (for Resolvent analysis)
        if gain_number is not None and self.analysisType == AnalysisType.RESOLVENT:
            matchingGainIndices         = [i for i, gNum in enumerate(gainNb) if gNum == gain_number]
            if len(matchingGainIndices) == 0:
                logger.error(f'No matching gain numbers found in the import folder for the specified gain number "{gain_number}".')
                return
            h5Files                     = [h5Files[i] for i in matchingGainIndices]
            logger.info(f'Importing only modes with gain number: {gain_number}')

        # Import each mode file
        numberOfModes                   = len(h5Files)
        logger.info(f'{numberOfModes} modes to import into mode collection.')
        for i in range(numberOfModes):
            logger.info(f'Importing mode {i+1}/{numberOfModes} from file "{h5Files[i]}".')
            
            # Get the mode type
            if self.analysisType == AnalysisType.MODAL:
                if 'Direct' in h5Files[i]:
                    modeSet = 'Direct'
                else:
                    modeSet = 'Adjoint'
            
            elif self.analysisType == AnalysisType.RESOLVENT:
                if 'Response' in h5Files[i]:
                    modeSet = 'Response'
                else:
                    modeSet = 'Forcing'
                    
            elif self.analysisType == AnalysisType.INPUT_OUTPUT:
                modeSet = 'Response'  # Always response for IO modes
            
            # TODO: find way to get the wavenumber "m", is always assumed to be zero otherwise
            mode                        = Mode(
                self.femSpace,
                self.mesh,
                isStateVector           = self.isStateVector,
                analysisType            = self.analysisType.name,
                modeType                = modeSet,
            )
            
            # Import the mode data from file
            mode.import_data(
                reader,
                importFolder,
                importFileName = h5Files[i],
            )
            logger.info(f'Imported mode {i+1}/{numberOfModes} from "{h5Files[i]}".')
            self.modeList.append(mode)
