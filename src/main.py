#/* FELICS, Finite Element Linearized Combustion Solver Copyright (C) 
# *2019 FLOW group TU Berlin - All Rights Reserved
# *
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally 
# * granted when asked for.
# * 
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
import 	argparse
from 	runModal       				import  runModal      
from 	runResolvent   				import  runResolvent  
from 	runInputOutput 				import  runInputOutput
from 	FELiCS.Parameters.config	import 	config
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
    """

# Initialize the argument parser
parser = argparse.ArgumentParser(description=desc_text, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
# Define the `-f` flag followed by a required argument for the config file
parser.add_argument('-f',"--file", type=str, required=True, metavar="path",help='Specify the path to the config file')
parser.add_argument('-d', '--debug', action='store_true', help='activate debug mode for extended output')
args = parser.parse_args()

# Initialize the logger
logger = Logger(args.debug, "felics")

if __name__ == '__main__':

    # Get the parameters
    param = config()
    param.importFromFile(args.file)

    mode = param.Case.AnalysisMode
    if mode == "Modal":
        runModal(param)
    elif mode == "Resolvent":	
        runResolvent(param)
    elif mode == "Input-Output":
        runInputOutput(param)