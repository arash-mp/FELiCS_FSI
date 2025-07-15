import  time
import  numpy                   as np
from 	FELiCS.Misc.logging     import Logger

# Get the logger
logger = Logger.get_logger("felics")

def runResolvent(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    from    FELiCS.IO.ExportSolution            import ExportFromFile 
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces
    from    FELiCS.Fields.meanFlowClass         import meanFlowClass
    from    FELiCS.Fields.fluctuationClass      import fluctuationSolutions
    from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
    from    FELiCS.Misc.functions               import printDebug
    from    FELiCS.Solvers.LinearSolver         import LinearSolver, ResolventOperator
    # import FELiCS.IO.Import as Import

    logger.info("Running Resolvent analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Get the mesh
    mesh = param.getMesh()
    
    # Define the FEM spaces
    FEMSpaces = FEMSpaces(
            param,
            mesh,
        )
        
    # read in mean flow and export to h5-file
    meanFlow = meanFlowClass(param, FEMSpaces, mesh)
    meanFlow.importDataFromFileAndExportToH5()
 
    # Define the equations
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

    P_shrink_f = equation.getShrinkerMatForcing()   
    P_forcing  = equation.getRestrictorMatForcing().matMult(P_shrink_f)

    P_shrink_r = equation.getShrinkerMatResponse()
    P_response = P_shrink_r.transposeMatMult(equation.getRestrictorMatResponse())

    W_f_big    = equation.getResolventNorm_forcing (meanFlow)
    W_r_big    = equation.getResolventNorm_response(meanFlow)

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

    # Start tracking time
    start = time.time()
    
    # TODO: use a Modecollection object here
    for i, omega in enumerate(omegas):

        logger.info("Solving resolvent SVD for omega = " + str(omega))

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
            X1, X2 = W_forcing.getVecs()             
            X1.setValues(range(0,len(eigenvectors_c[:,k])),eigenvectors_c[:,k])             
            Y1, Y2 = W_FEM.getVecs()             

            # forcings = Pu*eigenVectors             
            P_forcing.mult(X1,Y1)             
            forcings[:,k,i] = Y1.getValues(range(0,Y1.getSize()))             

            # Y1 = -1j * B_femWeight * forcings             
            W_FEM.mult(Y1,Y2)             
            Y2.scale(-1j)             

            # solve (A-omega*B)*responses = Y1             
            resolventOperator.getKSP().solve(Y2,Y1)             
            responses[:, k, i] = Y1.getValues(range(0, Y1.getSize()))


    # end tracking time
    end = time.time() - start
    logger.info('Solving the SVD(s) took %4g s' % end)

    #-----------------------------------------------------------------------
    ## EXPORT SOLUTION
    #-----------------------------------------------------------------------
    fluctSolutList = []     
    # construct a fluctuationSolutions Object for the Response, Forcing of     
    # every Frequency     
    # TODO: should be replaced by a more elegant command, preferably in the Modecollection object
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



    ExportFromFile(param,FEMSpaces,fluctSolutList,meanFlow)
