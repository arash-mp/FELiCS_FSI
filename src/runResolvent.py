import pdb
import numpy as np
import copy
import time

def runResolvent(param, useGUI):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
        useGUI: Boolean, True if program is run using GUI, False if run from
        terminal directly
    '''
    import FELiCS.IO.Import as Import
    from   FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile

    import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
    from   FELiCS.Fields.meanFlowClass import meanFlowClass
    from   FELiCS.Fields.fluctuationClass import fluctuationSolutions
    from   FELiCS.Equation.EquationCollection import EquationCollectionClass
    from   FELiCS.Misc.functions import printDebug

    from   FELiCS.Solvers.LinearSolver import LinearSolver, ResolventOperator 


    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh=param.BCs.getMesh()
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
    # read in mean flow
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFile()
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


    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # get matrices for resolvent analysis
    A          = equation.getLinearOperator(meanFlow)
    B          = equation.getWeightMatrix  (meanFlow)
    W_FEM      = equation.getResolventWeighting_FEM(meanFlow)  #=> including tensorial stuff

    P_forcing  = equation.getPMat()
    P_shrink_f = equation.getSimplePMat()   #shrinker mat for forcing=> write routine

    P_shrink_r = equation.getShrinkerMat_response()
    P_response = P_shrink_r.transposeMatMult(equation.getCrMat())

    W_f_big    = equation.getResolventNorm_forcing (meanFlow)
    W_r_big    = equation.getResolventNorm_response(meanFlow)

    #W_forcing  = P_shrink_f.transposeMatMult( W_FEM.matMult(P_shrink_f))
    #W_forcing  = P_forcing.transposeMatMult(W_FEM.matMult(P_forcing))

    W_forcing  = P_shrink_f.transposeMatMult(W_f_big.matMult(P_shrink_f))
    W_response = P_shrink_r.transposeMatMult(W_r_big.matMult(P_shrink_r))


    # get parameters for resolvent analysis
    omegas   = param.IOResolvent.Omegas
    nSol     = param.Numerics.nSolut
    nOmegas  = len(omegas)
    nDofs    = A.getSizes()[0][0]


    gains     = np.zeros((nSol, nOmegas),'complex')
    forcings  = np.zeros((nDofs, nSol, nOmegas),'complex')
    responses = np.zeros((nDofs, nSol, nOmegas),'complex') 


    # track time
    start= time.time()

    for i, omega in enumerate(omegas):

        printDebug(True, "-- Performing resolvent analysis for omega = " + str(omega))

        # R = A-omega*B         
        R = A.copy()         
        R.axpy(-omega, B)            
        resolventOperator = ResolventOperator(                                         
                R,                                         
                W_FEM,                                         
                W_forcing,                                         
                W_response,                                         
                P_forcing,                                         
                P_response)
        # Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"         
        # via the matrix vector multiplation "mult"
        gains[:,i],eigenvectors_c = LinearSolver.solveSVDOfResolvent(
                resolventOperator,
                nev=nSol,
                tol=1.e-16,
                max_it=200,
                )
        # Write gains to results dictionary
        gains[:, i] = np.real(gains[:, i])

        # Iterate through the first nSolut gains
        # Compute the respetive forcing and responses with the solution of the SVD ("eigenvetors_c")
        # TODO Sophie: do this more elegantly
        for k in range(nSol):
            
            # get petsc vectors from petsc matrices (=get petsc vectors with correct sizes)             
            X1, X2 = P_forcing.getVecs()             
            X1.setValues(range(0,len(eigenvectors_c[:,k])),eigenvectors_c[:,k])             
            Y1, Y2 = W_FEM.getVecs()             

            # forcings = Pu*eigenVectors             
            P_forcing.mult(X1,X2)             
            forcings[:,k,i] = X2.getValues(range(0,X2.getSize()))             

            # Y1 = -1j * B_femWeight * forcings             
            W_FEM.mult(X2,Y1)             

            # Y1 = -1j * B_femWeight * forcings             
            W_FEM.mult(X2,Y1)             
            Y1.scale(-1j)             

            # solve (A-omega*B)*responses = Y1             
            resolventOperator.getKSP().solve(Y1,X2)             
            responses[:, k, i] = X2.getValues(range(0, X2.getSize()))


    # end tracking time
    end = time.time() - start
    printDebug(True, '-- Solving the SVD(s) for the resolvent took %4g s' % end)
    #residuum_max = solution.getMaximumError()
    #printDebug(True, '-- Maximum residuum of all solutions:  %12g' % (residuum_max))



    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    #fluctSolutList = solution.getOldSolutionObject(meanFlow, param, FEMSpaces)

    fluctSolutList = []     
    # construct a fluctuationSolutions Object for the Response, Forcing of     
    # every Frequency     
    for i, omega in enumerate(omegas):         
        for gainNumb in range(responses.shape[1]):             
            fluctSolutForcing = fluctuationSolutions(                                 
                    param,                                 
                    meanFlow,                                 
                    FEMSpaces,                                 
                    omegas[i],                                 
                    forcings[:,gainNumb,i],                                 
                    False,                                 
                    gainNumb,                                 
                    gains[gainNumb, i],                                 
                    )             

            fluctSolutResponse = fluctuationSolutions(                                 
                    param,                                 
                    meanFlow,                                 
                    FEMSpaces,                                 
                    omegas[i],                                 
                    responses[:,gainNumb, i],                                 
                    True,                                 
                    gainNumb,                                 
                    gains[gainNumb, i],                                 
                    )             
            fluctSolutList.append(fluctSolutForcing)             
            fluctSolutList.append(fluctSolutResponse)



    if useGUI:
        ExportGUI(param, fluctSolutList, meanFlow,FEMSpaces, equation,mesh)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)

