import pdb
import numpy as np
import copy
import time
import sys
import dolfinx

from   FELiCS.Parameters.parameters import parameters

import FELiCS.IO.Import as Import
import FELiCS.IO.ExportSolution as Export 

import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
from   FELiCS.Fields.meanFlowClass import meanFlowClass
from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
from   FELiCS.Equation.EquationCollection import EquationCollectionClass
from   FELiCS.Misc.functions import printDebug

from   FELiCS.Solvers.LinearSolver import LinearSolver 
from   FELiCS.Fields.ModeCollection import ModeCollection
from   FELiCS.Fields.Field import Field
from   FELiCS.Fields.Mode import Mode
from   FELiCS.Misc.tensorUtils import Tensor

from   CaseHandler import CaseHandler


def calculateResolventModesAndSensitivity(settingsFileName, baseFlow_array, baseFlowSensitivity, optimizerParameters, deformed = False):

    #-----------------------------------------------------------------------
    ## INITIALIZATION 
    #-----------------------------------------------------------------------
    # read parameters
    param=parameters()
    param.importFromFile(settingsFileName)
    param.getOldParameters()

    # mesh
    mesh=param.BCs.getMesh()
   
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
    
    
    # initialize mean flow class
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    
    #import mean flow data from file
    meanFlow.importDataFromFile()

    # load base flow into meanFlow object
    if np.linalg.norm(meanFlow._fieldDict['u'].x.array[:]) < 1.e-8:
        baseFlow = Field(FEMSpaces.VMixed, mesh)
        baseFlow.setCoefficientArray(baseFlow_array)
        [u,p] = baseFlow.getListOfSingleFields()
        meanFlow._fieldDict['u'] = u.function

    # export mean flow in "h5" file
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        meanFlow.exportBaseFlowAsHDF5()
    meanflowFilename = 'meanflow.h5'
    meanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)

    # equation
    equation = EquationCollectionClass(
                                      param,
                                      FEMSpaces,
                                      meanFlow,
                                      mesh
                                      )

    caseHandler = CaseHandler(settingsFileName, param, mesh, equation.boundaries)


    # deform mesh at beginning - only for testing /debugging
    if deformed == True:
        geometryDeformer = caseHandler.getGeometryDeformer()
        geometryDeformer.deformMesh(optimizerParameters)
        geometryDeformer.isDeformed = False


    #-----------------------------------------------------------------------
    ## SOLVE RESOLVENT EIGENPROBLEM AND SCALE MODES 
    #-----------------------------------------------------------------------
    # get matrices for eigenproblem
    A  = equation.getLinearOperator(meanFlow)
    B  = equation.getWeightMatrix  (meanFlow)
    
    # get parameters for eigenproblem
    omegas   = param.IOResolvent.Omegas
    nSol     = param.Numerics.nSolut
    
    # track time
    start= time.time()
    
    printDebug(True, "-------------------------------------------------------------" )
    printDebug(True, "-- Calculating dominant resolvent modes and gains...")
    printDebug(True, "-------------------------------------------------------------" )


#################################################################################################
#################################################################################################
################# "CONSTRUCTION SITE" ###########################################################
#################################################################################################
#################################################################################################
    Restrictor_forcing  = equation.getPMat()
    Restrictur_response = equation.getCrMat()

    WeightMatrix        = equation.getFEMWeightMatrix()

    exit()

    self.__matrix_dict_petsc['Q'] = \
            self.__matrix_dict_petsc['Pu'].transposeMatMult(\
                    WeightMatrix.matMult(\
                    self.__matrix_dict_petsc['Pu']))

    for omega in omegas:

        printDebug(True, "-- Performing resolvent analysis for omega = " + str(omega))

        OP_petsc          = self.__matrix_dict_petsc['A'].copy()
        OP_petsc.axpy(-omega, self.__matrix_dict_petsc['B'])    # OP_petsc = A-omega*B

        resolventOperator = ResolventOperator(
                                        OP_petsc, 
                                        self.__matrix_dict_petsc['Q'], 
                                        self.__matrix_dict_petsc['Pu'],
                                        self.__matrix_dict_petsc['Cr'],
                                        self.__matrix_dict_petsc['B_forcing'],
                                        self.__matrix_dict_petsc['B_response'])
    
        # Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"
        # via the matrix vector multiplation "mult"
        gains[:,i],eigenvectors_c = LinearSolver.solveSVDOfResolvent(
                            resolventOperator,
                            nev=nSol,
                            tol=1.e-13,
                            max_it=200,
                            )

        # Write gains to results dictionary
        gains[:, i] = np.real(gains[:, i])

        # Iterate through the first nSolut gains
        # Compute the respetive forcing and responses with the solution of the SVD ("eigenvetors_c")
                        # TODO Sophie: do this more elegantly
        for k in range(self.__param.Numerics.nSolut):

            # get petsc vectors from petsc matrices (=get petsc vectors with correct sizes)
            X1, X2 = self.__matrix_dict_petsc['Pu'].getVecs()
            X1.setValues(range(0,len(eigenvectors_c[:,k])),eigenvectors_c[:,k])
            Y1, Y2 = self.__matrix_dict_petsc['B_femWeight'].getVecs()

            # forcings = Pu*eigenVectors
            self.__matrix_dict_petsc['Pu'].mult(X1,X2)
            forcings[:,k,i] = X2.getValues(range(0,X2.getSize()))

            # Y1 = -1j * B_femWeight * forcings
            self.__matrix_dict_petsc['B_femWeight'].mult(X2,Y1)

            # Y1 = -1j * B_femWeight * forcings
            self.__matrix_dict_petsc['B_femWeight'].mult(X2,Y1)
            Y1.scale(-1j)

            # solve (A-omega*B)*responses = Y1
            resolventOperator.getKSP().solve(Y1,X2)
            responses[:, k, i] = X2.getValues(range(0, X2.getSize()))

            resolventOperator.destroySelf()

        # construct a fluctuationSolutions Object for the Response, Forcing of
        # every Frequency
        for i, omega in enumerate(self.__param.IOResolvent.Omegas):
            for gainNumb in range(responses.shape[1]):
                fluctSolutForcing = fluctuationSolutions(
                                    self.__param,
                                    self.__meanFlow,
                                    self.__FEMSpaces,
                                    self.__param.IOResolvent.Omegas[i],
                                    forcings[:,gainNumb,i],
                                    False,
                                    gainNumb,
                                    gains[gainNumb, i],
                                    )
                fluctSolutResponse = fluctuationSolutions(
                                    self.__param,
                                    self.__meanFlow,
                                    self.__FEMSpaces,
                                    self.__param.IOResolvent.Omegas[i],
                                    responses[:,gainNumb, i],
                                    True,
                                    gainNumb,
                                    gains[gainNumb, i],
                                    )

                fluctSolutObjList.append(fluctSolutForcing)
                fluctSolutObjList.append(fluctSolutResponse)

        toc_res = time.perf_counter() - tic_res
        printDebug(True, f"-- Solving resolvent took: {toc_res:0.4f} seconds")

        return fluctSolutObjList


#################################################################################################
#################################################################################################


    # get leading eigenvalue
    eigenValue    = solution.getLeadingMode().getEigenValue()
     
    # solve adjoint eigenproblem only for the leading eigenvalue
    tmp = LinearSolver.solveGeneralEigenproblem(A,
                                                B,
                                                eigenValue,
                                                1,
                                                adjoint=True)
    solution.appendSolutionOfEigenProblem(tmp, guess, adjoint=True)
    
    
    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the resolvent eigenproblem took %4g s' % end)
    residuum_max = solution.getMaximumError()
    printDebug(True, '-- Maximum residuum of all solutions:  %12g' % (residuum_max))
    
    # get leading modes
    mode_direct  = solution.getLeadingMode(adjoint=False)
    mode_adjoint = solution.getLeadingMode(adjoint=True)
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- Leading eigenvalue:  ' + str(eigenValue))
    printDebug(True, '------------------------------------------------ ')
    
    # scale modes s.t. mode_adjoint^H * B * mode_direct = 1
    mode_direct_petsc  = mode_direct.getPetscVector()
    mode_adjoint_petsc = mode_adjoint.getPetscVector()
    temp               = mode_direct.getPetscVector() # gets a petsc vector "temp" of correct length
    B.mult(mode_direct_petsc, temp)                   # temp = B*mode_direct
    factor = temp.dot(mode_adjoint_petsc)             # factor = mode_adjoint ^H temp
    
    mode_direct_petsc.scale(1./np.sqrt(factor))
    mode_adjoint_petsc.scale(np.conj(1./np.sqrt(factor)))
    mode_direct.setCoefficientArray(mode_direct_petsc.getArray())
    mode_adjoint.setCoefficientArray(mode_adjoint_petsc.getArray())

    ## export all modes and all eigenvalues in standard felics format
    #fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)
    #ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

    # export leading modes in standard felics format
    solution_onlyLeading = ModeCollection(FEMSpaces.VMixed, mesh)
    solution_onlyLeading.appendMode(mode_direct)
    solution_onlyLeading.appendMode(mode_adjoint)
    fluctSolutList = solution_onlyLeading.getOldSolutionObject(meanFlow, param, FEMSpaces)
    Export.ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

    # export full (direct) spectrum (and overwrite the spectrum from the other export function, which only writes the leading eigenvalue)
    spectrum_direct = solution.getDirectEigenValueSpectrum()
    Export.writeCSVSpectrum(param,spectrum_direct)


    #-----------------------------------------------------------------------
    ## CALCULATE GRADIENT WITH RESPECT TO GEOMETRY DEFORMATIONS - Part I
    #-----------------------------------------------------------------------
    # get undisturbed operators (#TODO: cannot use the ones from above - why?)
    A_0  = equation.getLinearOperator(meanFlow)
    B_0  = equation.getWeightMatrix  (meanFlow)

    # parameters that are used in the optimization process
    a_i  = optimizerParameters 
    
    # epsilon for mesh deformation
    epsilon = 1.e-8
   
    # initialize stuff     
    geometryDeformer = caseHandler.getGeometryDeformer() 
    N_param          = geometryDeformer.getNumberOfParameters()
    sensitivity1     = np.zeros(N_param,dtype=complex)
    sensitivity2     = np.zeros(N_param,dtype=complex)
    baseFlowSens     = Field(FEMSpaces.VMixed, mesh)

    for i in range(0,N_param):
        printDebug(True, "-------------------------------------------------------------" )
        printDebug(True, "-- Calculating eigenvalue sensitivity to parameter number "+str(i+1)+"...")
        printDebug(True, "-------------------------------------------------------------" )

        a_i[i] = a_i[i] + epsilon
    
        geometryDeformer.deformMesh(a_i)
       
        ## calculate the first sensitivity part: with the partial derivative of the linear operator and the weight matrix
        # (A_deformed - A_0)/epsilon
        A_deformed  = equation.getLinearOperator(meanFlow)
        A_deformed.axpy(-1., A_0)    # A_deformed -= A_0
        A_deformed.scale(1./epsilon) # A_deformed /= epsilon

        # (B_deformed - B_0)/epsilon 
        B_deformed  = equation.getWeightMatrix  (meanFlow)
        B_deformed.axpy(-1., B_0)     # B_deformed -= B_0
        B_deformed.scale(1./epsilon)  # B_deformed /= epsilon 

        # multiply with eigenvalue
        A_deformed.axpy(-eigenValue, B_deformed) # A_deformed -= eigenValue * B_deformed

        # multiply with the adjoint eigenvector (from left) and the direct eigenvector (from right)  
        mode_direct_petsc  = mode_direct.getPetscVector()
        mode_adjoint_petsc = mode_adjoint.getPetscVector()
        result             = mode_direct.getPetscVector()
        A_deformed.mult(mode_direct_petsc, result)
        dolfinx.fem.petsc.set_bc(result, equation.BCs)
        sensitivity1[i]   = result.dot(mode_adjoint_petsc)
       
        ## calculate the second sensitivity part: with the base flow sensitivities
        baseFlowSens.setCoefficientArray(baseFlowSensitivity[i])
        [u_sens, p_sens]                  = baseFlowSens.getListOfSingleFields()
        meanFlow._fieldDict['u_bilinear'] = u_sens.function
        BL                                = equation.getBilinearOperator(meanFlow)
        BL.mult(mode_direct_petsc, result)
        #dolfinx.fem.petsc.set_bc(result, equation.BCs)
        sensitivity2[i]  = result.dot(mode_adjoint_petsc)


        a_i[i] = a_i[i] - epsilon
        geometryDeformer.restoreMesh()
    
    
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- sensitivities part 1:  ' + str((sensitivity1)))
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '------------------------------------------------ ')
    printDebug(True, '-- sensitivities part 2:  ' + str((sensitivity2)))
    printDebug(True, '------------------------------------------------ ')


    return sensitivity1, sensitivity2, eigenValue 


    



class ResolventOperator(object):

        """
        This class serves as a "matrix-free" representation of the Resolvent operator multiplicated with its Hermitian transposed, to conduct the singular value decomposition of the system. 
        It contains a method called "mult", which is called by the eigenvalue solver, and returns a matrix vector product of the represented matrix. 

        Private attributes:
        Protected attributes:
        Public attributes:
        """

        def __init__(self,
                OP,
                Qf,
                Pu,
                Cr,
                B_forcing,
                B_response):


                from petsc4py import PETSc

                self._size = Qf.getSize()

                self._Z1, self._Z2  = OP.getVecs()
                self._Y1, self._Y2  = Qf.getVecs()

                self._Pu = Pu
                self._Cr = Cr
                self._Bf = B_forcing
                self._Br = B_response

                # create KSP1: This is a solver for the System OP*x=y.
                self._ksp1 = PETSc.KSP().create()
                self._ksp1.setOperators(OP)
                self._ksp1.setType(PETSc.KSP.Type.PREONLY)
                self._ksp1.getPC().setType(PETSc.PC.Type.LU)
                self._ksp1.getPC().setFactorSolverType('mumps')
                self._ksp1.setUp()
                # create KSP2: This is a solver for the System conj(OP)*x=y. It will later be used to solve the transposed system, thus effectively solving OP^H *x=y, which is the Hermitian transpose of the system.
                # TODO Sophie: unfortuantely there is no "solveHermitianTranspose" in the petsc4py (yet?). Thus we have to do an additional LU decomposistion.... Change as soon as this is included in the petsc4py! 
                OP_H = OP.copy()   #create a new matrix, s.t. the original one will not be overwritten 
                OP_H.conjugate()
                OP_H.assemble()
                self._ksp2 = PETSc.KSP().create()
                self._ksp2.setOperators(OP_H)
                self._ksp2.setType(PETSc.KSP.Type.PREONLY)
                self._ksp2.getPC().setType(PETSc.PC.Type.LU)
                self._ksp2.getPC().setFactorSolverType('mumps')
                self._ksp2.setUp()
                # create KSP3: This is a solver for the System Qf*x=y (will later be used to solve the transposed system).
                #Qf.conjugate() #=> is this needed?
                self._ksp3 = PETSc.KSP().create()
                self._ksp3.setOperators(Qf)
                self._ksp3.setType(PETSc.KSP.Type.PREONLY)
                self._ksp3.getPC().setType(PETSc.PC.Type.LU)
                self._ksp3.getPC().setFactorSolverType('mumps')
                self._ksp3.setUp()

        def getSize(self):
                return self._size

        def getVecs(self):
                return self._Y1, self._Y2

        def mult(self, mat, X, Y):
                # returns Y=mat*X 
                # mat = (Qf^T)^-1 * Pu^T * Bf^T * (OP^H)^-1 * Cr^T * Br * Cr * OP^-1 * Bf * Pu

                self._Pu.mult            (X,        self._Z1)  #Z1 = Pu*X
                self._Bf.mult            (self._Z1, self._Z2)  #Z2 = Bf*Z1
                self._ksp1.solve         (self._Z2, self._Z1)  #Z1 = OP^-1 * Z2
    
                self._Cr.mult            (self._Z1, self._Z2)  #Z2 = Cr*Z1
                self._Br.mult            (self._Z2, self._Z1)  #Z1 = Br*Z2
                self._Cr.multTranspose   (self._Z1, self._Z2)  #Z2 = Cr^T * Z1
                self._ksp2.solveTranspose(self._Z2, self._Z1)  #Z1 = (OP^H)^-1 * Z2

                self._Bf.multTranspose   (self._Z1, self._Z2)  #Z2 = Bf^T * Z1
                self._Pu.multTranspose   (self._Z2, self._Y1)  #Y1 = Pu^T * Z2
                self._ksp3.solveTranspose(self._Y1, Y)         #Y  = (Qf^T)^-1 * Y1

                return Y

        def getKSP(self):
                return self._ksp1

        def destroySelf(self):
                """Clean-up the disc space to avoid memory leaks. Should be called if several resolvent SVDs are done one after the other.""" 
                self._ksp1.getPC().destroy()
                self._ksp2.getPC().destroy()
                self._ksp3.getPC().destroy()
                self._ksp1.destroy()
                self._ksp2.destroy()
                self._ksp3.destroy()




