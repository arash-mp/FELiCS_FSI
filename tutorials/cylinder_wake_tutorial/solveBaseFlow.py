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
import pdb
import numpy as np
import copy
import time
import h5py

import dolfinx

from   mpi4py import MPI

from   FELiCS.Parameters.Config import Config

from   FELiCS.IO.Writer                   import Writer
from   FELiCS.SpaceDisc.FEMSpaces         import FEMSpaces
from   FELiCS.Fields.MeanFlowClass        import MeanFlowClass
from   FELiCS.Equation.EquationCollection import EquationCollectionClass

from   FELiCS.Solvers.LinearSolver  import LinearSolver 
from   FELiCS.Fields.Field          import Field


def createInitialSolutionAsFelFile():
    xx = np.linspace(-101, 201, 400)
    yy = np.linspace(-26, 26, 100)
    X, Y = np.meshgrid(xx, yy, indexing='ij')
    x = (X.flatten())
    y = (Y.flatten())
    
    r = np.sqrt(x**2+y**2)
    ux = np.zeros_like(x)
    ux[r>2] = 1 - np.exp(-((r[r>2]-1)**2)/0.1)
    ux[r<=2] = 0
    
    file = "initial_solution.fel"
    
    with h5py.File(file, 'w') as f:
        f.create_dataset('/MeanFlow/x',   data=x)
        f.create_dataset('/MeanFlow/y',   data=y)
        f.create_dataset('/MeanFlow/ux',  data = ux)


def createFelFileFromH5(felFileName, h5FileName, meshFileName):
    with h5py.File(h5FileName, 'r') as f:
        ux = f['ux'][()]
        uy = f['uy'][()]
        p  = f['p'][()]
    with h5py.File(meshFileName, 'r') as f:
        x = f['/coordinates/x'][()]
        y = f['/coordinates/y'][()]
    with h5py.File(felFileName, 'w') as f:
        f.create_dataset('/MeanFlow/x',   data=x)
        f.create_dataset('/MeanFlow/y',   data=y)
        f.create_dataset('/MeanFlow/ux',  data = ux)
        f.create_dataset('/MeanFlow/uy',  data = uy)
        f.create_dataset('/MeanFlow/p',   data = p)


def calculateBaseFlow(settingsFileName, optimizerParameters = None, deformed = False):

    #-----------------------------------------------------------------------
    ## INITIALIZATION 
    #-----------------------------------------------------------------------
    # read parameters
    param=Config()
    param.import_from_file(settingsFileName)

    # mesh
    mesh      = param.get_mesh()
    
    # FEMSpaces
    femSpaces = FEMSpaces(param, mesh)
    
    # writer to export the results in files
    writer    = Writer(mesh, param.Export.ExportFolder)

    # initialize mean flow class 
    meanFlow = MeanFlowClass(param, femSpaces, mesh)
    meanFlow.import_data_from_file_and_export_to_h5(writer)
    
    # equation
    equation = EquationCollectionClass(
                                      param,
                                      femSpaces,
                                      meanFlow,
                                      mesh
                                      )

    #-----------------------------------------------------------------------
    ## CALCULATE BASE FLOW 
    #-----------------------------------------------------------------------
    
    # Define a Field class to save calculating results
    baseFlow     = Field(femSpaces.VMixed, mesh, isStateVector=True)
    # Initialize fields in baseFlow
    mapping_ux_b = baseFlow.space.sub(0).sub(0).collapse()[1]
    mapping_ux_m = meanFlow._fieldDict["u"].space.sub(0).collapse()[1]
    baseFlow.function.x.array[mapping_ux_b] = meanFlow._fieldDict["u"].function.x.array[mapping_ux_m]

    # ---------------- START LOOP ---------------------------------------------- 
    ## start Newton solver
    target_residuum = 1.e-11
    # track time
    start= time.time()
    # calculate initial residuum
    N = equation.get_nonlinear_expression(meanFlow)
    residuum = np.linalg.norm(N.getArray())
    # calulate the base flow
    i=0
    while(residuum > target_residuum and i<100):
        i+=1

        # solve equation system
        L = equation.get_linear_operator(meanFlow)
        newtonSummand_array = LinearSolver.solve_equation_system(L,N)
        L.destroy()
        N.destroy()

        # update baseFlow
        baseFlow_array = baseFlow.get_coefficient_array() - newtonSummand_array
        baseFlow.set_coefficient_array(baseFlow_array)

        # calculate nonlinear expression & residuum
        [u,p] = baseFlow.get_list_of_sub_fields()
        meanFlow._fieldDict['u'] = u
        meanFlow._fieldDict['p'] = p
        N = equation.get_nonlinear_expression(meanFlow)
        residuum = np.linalg.norm(N.getArray())
        print("-------------------------------------------------------------" )
        print("-- Base flow iteration: "+str(i)+"; Residuum: %4g " % residuum)
        print("-------------------------------------------------------------" )

    
    # end tracking time
    end = time.time() - start
    print('-- Solving the base flow problem took %4g s' % end)
    print('-- Residuum:  %12g' % (residuum))

    baseFlow.export_to_h5(writer, "baseFlow")
 


##################################################################################################
# MAIN SRIPT

# 1. create an initial solution for the base flow computation
createInitialSolutionAsFelFile()

# 2. calculate the base flow
configFilename       = 'Re50.json'
calculateBaseFlow(configFilename)

# 3. create a fel file from the result, which can be used as input again
createFelFileFromH5("base_flow_for_FELiCS.fel", "out/baseFlow.h5", "out/mesh.h5")
