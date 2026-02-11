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
import 	argparse
from 	run_modal     				import  run_modal  
from 	run_resolvent    			import  run_resolvent
from 	run_input_output  			import  run_input_output
# from 	runModal_testReader     	import  runModal      
# from 	runResolvent_testReader    	import  runResolvent  
# from 	runInputOutput_testReader  	import  runInputOutput
from 	FELiCS.Parameters.Config	import 	Config
from 	FELiCS.Misc.logging			import  Logger

# Define the description and epilog for the help message
desc_text = """
------------------ FELICS -----------------
Finite Element Linearized Combustion Solver
**add short FELICS description here**
"""
epilog = """
Example usage:
    python main.py -h		shows help message 
    python main.py -f path	start with config file at path (mandatory)
    python main.py -f path -d	for debug mode
    python main.py -f path -t	for test mode
    """

# Initialize the argument parser
parser = argparse.ArgumentParser(
description=desc_text,
epilog=epilog,
formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
'-f',
"--file",
"-file",
type=str,
required=True,
metavar="path",
help='Specify the path to the config file',
)
parser.add_argument(
'-d',
'--debug',
action='store_true',
help='activate debug mode for extended output',
)
parser.add_argument(
'-t',
'--test',
action='store_true',
help='activate test mode with no output',
)
args = parser.parse_args()

# Initialize the logger
logger = Logger(
args.debug,
args.test,
"felics",
)
logger = Logger.get_logger("felics")

if __name__ == '__main__':

    # Get the parameters
    param = Config()
    param.import_from_file(args.file)
    mode = param.Case.AnalysisMode
    if mode == "Modal":
        run_modal(param)
    elif mode == "Resolvent":
        run_resolvent(param)
    elif mode == "Input-Output":
        run_input_output(param)
    else:
        logger.error("Analysis type not recognized in FELiCS main.")
    logger.info("Finished FELiCS run.")
