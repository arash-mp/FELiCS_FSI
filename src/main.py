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


"""
This file is called to start the program

This file was created by Thomas L. Kaiser. Significant contributions 

Parameters
----------


""" 

import tkinter

import sys

from runModal       import  runModal      
from runResolvent   import  runResolvent  
from runInputOutput import  runInputOutput

from FELiCS.GUI.GUI import FELiCS_GUI
from FELiCS.Parameters.parametersClass import parameters

# check if '-file' argument was added to run from file only. else start the GUI

# Replace spurious double quotation marks which occur in MobaXterm when executing shell scripts
for i_argument,argument in enumerate(sys.argv):
        sys.argv[i_argument] = argument.replace(chr(8221), '')

if len(sys.argv) == 1:
    useGUI= True
elif sys.argv[1] == '-file': useGUI=False
else: useGUI=True

if __name__ == '__main__':
    if useGUI: # Run program in GUI mode
        window=FELiCS_GUI()
    else: # Run program in terminal mode from settings file
        SettingsFileName = sys.argv[2]
        param=parameters()
        param.importFromFile(SettingsFileName)
        param.getOldParameters()

        mode = param.Case.AnalysisMode
        if mode == "Modal":
                runModal(param,useGUI=False)
        elif mode == "Resolvent":
                runResolvent(param,useGUI=False)
        elif mode == "Input-Output":
                runInputOutput(param,useGUI=False)
