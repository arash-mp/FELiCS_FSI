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
import  time

# Local Libraries and methods
from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
from    FELiCS.Fields.MeanFlowClass         import MeanFlowClass
from    FELiCS.Fields.ModeCollection        import ModeCollection
from    FELiCS.IO.Writer                    import Writer
from 	FELiCS.Misc.logging                 import Logger
from    FELiCS.Solvers.LinearSolver         import LinearSolver, ResolventOperator

# Get the logger
logger                          = Logger.get_logger("felics")

def run_resolvent(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    logger.warning("Running Resolvent analysis")
    
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # mesh
    mesh      = param.get_mesh()
    from    FELiCS.SpaceDisc.FEMSpaces          import FEMSpaces

    # FEMSpaces
    FEMSpaces = FEMSpaces(
        param,
        mesh,
    )
     
    # writer to export the results in files
    writer    = Writer(
        mesh,
        param.Export.ExportFolder,
    )

    # read in mean flow and export to h5-file
    meanFlow = MeanFlowClass(
        param,
        FEMSpaces,
        mesh,
    )
    meanFlow.import_data_from_file_and_export_to_h5(writer)

    # equation
    equation = EquationCollectionClass(
        param,
        FEMSpaces,
        meanFlow,
        mesh,
    )

    #-----------------------------------------------------------------------
    ## MAIN PART
    #-----------------------------------------------------------------------
    # Get operators required for resolvent analysis
    A                           = equation.get_linear_operator(meanFlow)
    B                           = equation.get_weight_matrix  (meanFlow)
    W_FEM                       = equation.get_resolvent_weighting_fem(meanFlow)  #=> including tensorial stuff

    P_shrink_f                  = equation.get_shrinker_mat_forcing()   
    P_forcing                   = equation.get_restrictor_mat_forcing().matMult(P_shrink_f)

    P_shrink_r                  = equation.get_shrinker_mat_response()
    P_response                  = P_shrink_r.transposeMatMult(equation.get_restrictor_mat_response())

    W_f_big                     = equation.get_resolvent_norm_forcing (meanFlow)
    W_r_big                     = equation.get_resolvent_norm_response(meanFlow)

    W_forcing                   = P_shrink_f.transposeMatMult(W_f_big.matMult(P_shrink_f))
    W_response                  = P_shrink_r.transposeMatMult(W_r_big.matMult(P_shrink_r))

    # get parameters for resolvent analysis
    omegas                      = param.IOResolvent.Omegas
    nSol                        = param.Numerics.nSolut
    
    # Start tracking time
    start                       = time.time()
    
    # Prepare solution object
    solution                    = ModeCollection(
        FEMSpaces.VMixed,
        mesh,
        analysisType='Resolvent',
    )
    
    # Run analysis at each frequency
    for i, omega in enumerate(omegas):

        logger.info("Solving resolvent SVD for omega = " + str(omega))
        
        # Catch bug for omega = 1.0
        if omega == 1.0:
            omega += 1.e-4
            logger.warning("Omega was equal to 1.0, which can lead to numerical issues. Added 1.e-4 to omega.")

        # initialize the resolvent operator, which is a class that imitates 
        # a matrix to use matrix-free methods
        R                       = A.copy()         
        R.axpy(
            -omega,
            B,
        )       # R = A-omega*B      
        resolventOperator       = ResolventOperator(
            R,
            W_FEM,
            W_forcing,
            W_response,
            P_forcing,
            P_response,
        )
        
        # Perform eigenvalue decomposition of the linear operator defined in the class "ResolventOperator"         
        # via the matrix vector multiplation "mult"
        gains, eigenvectors_c   = LinearSolver.solve_svd_of_resolvent(
            resolventOperator,
            nev             = nSol,
            tol             = 1.e-16,
            max_it          = 200,
        )

        solution.append_solution_of_svd_problem(
            eigenvectors_c,
            omega,
            gains,
            resolventOperator,
        )

        # export spectrum and newly calculated modes
        solution.export_spectrum_to_csv(writer)
        solution.export_modes(
            writer,
            onlyNewN = nSol*2,
        )
        
    # End tracking time
    logger.info(f"Solving the SVD(s) took {time.time()-start:.4g} s")


