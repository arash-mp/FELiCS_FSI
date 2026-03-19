# Standard libraries
import  time

# Third party libraries
import  numpy as np

# Local Libraries and methods
from    FELiCS.Equation.EquationCollection  import EquationCollectionClass
from    FELiCS.Fields.MeanFlowClass         import MeanFlowClass
from    FELiCS.Fields.ModeCollection        import ModeCollection
# import  FELiCS.IO.Import as Import
from    FELiCS.IO.Writer                    import Writer
from 	FELiCS.Misc.logging                 import Logger
from    FELiCS.Solvers.LinearSolver         import LinearSolver 

# Get the logger
logger = Logger.get_logger("felics")

def run_modal(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    logger.info("Running Modal analysis")
    #-----------------------------------------------------------------------
    ## INITIALIZATION
    #-----------------------------------------------------------------------
    # Get FELiCS objects required for analysis
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
    # get matrices for eigenproblem
    A = equation.get_linear_operator(meanFlow)
    B = equation.get_weight_matrix  (meanFlow)

    # get parameters for eigenproblem
    guesses          = param.Numerics.EigenValueGuess
    nSol             = param.Numerics.nSolut
    adjoint          = param.Case.CalculateAdjoint

    # track time
    start= time.time()

    # solve eigenproblem for each guess
    solution = ModeCollection(
        FEMSpaces.VMixed,
        mesh,
        analysisType = "modal",
    )
    for guess in guesses:
        
        logger.info("Solving direct GEVP for guess: omega = " + str(guess))
        tmp     = LinearSolver.solve_general_eigenproblem(
            A,
            B,
            guess,
            nSol,
        )

        solution.append_solution_of_eigen_problem(
            tmp,
            guess,
        )

        n_calculated = nSol

        if adjoint:
            
            logger.info("Solving adjoint GEVP for guess: omega = " + str(np.conj(guess)))
            tmp = LinearSolver.solve_general_eigenproblem(
                A,
                B,
                guess,
                nSol,
                adjoint=True,
            )

            solution.append_solution_of_eigen_problem(
                tmp,
                guess,
                adjoint=True,
            )

            n_calculated += nSol

        # Exporting the (temporary) spectrum to a file
        solution.export_spectrum_to_csv(writer)

        # Exporting only the newly calculated modes to files (for this guess)
        solution.export_modes(
            writer,
            onlyNewN = n_calculated,
        )


    # end tracking time
    end = time.time() - start
    logger.info('Solving the general eigenproblem took %4g s' % end)
    residuum_max    = solution.get_maximum_error()
    logger.debug('Maximum residuum of all solutions:  %12g' % (residuum_max))
   

