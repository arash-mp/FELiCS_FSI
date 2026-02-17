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
from    FELiCS.Solvers.LinearSolver         import LinearSolver 

# Get the logger
logger                      = Logger.get_logger("felics")

def run_input_output(param):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
    '''

    logger.info("Running input/output analysis")
    
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
    # Get operators required for input/output analysis
    A                       = equation.get_linear_operator(meanFlow)
    B                       = equation.get_weight_matrix  (meanFlow)
    forcing                 = equation.get_forcing_for_input_output(meanFlow) 

    # Get parameters for input/output analysis
    omegas                  = param.IOResolvent.Omegas

    # Start tracking time
    start                   = time.time()

    # Solve equation system for each frequency 
    solution                = ModeCollection(
    FEMSpaces.VMixed,
    mesh,
    analysisType = "input_output",
    )
    for omega in omegas:
        # define operator
        operator            = A.copy()
        operator.axpy(
        -omega,
        B,
        ) #petsc command: operator = A - omega*B
        solutionVector      = LinearSolver.solve_equation_system(
        operator,
        forcing,
        )
        solution.append_mode_from_vector(
        solutionVector,
        frequency = omega,
        gain = 1,
        ) 

        # export newest mode
        # TODO: there is no need for a ModeCollection, export modes directly
        solution.export_modes(
        writer,
        onlyNewN = 1,
        )

    # End tracking time
    logger.info(f"Solving the input/output problem took {time.time()-start:.4g} s")

