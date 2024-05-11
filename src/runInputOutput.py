import pdb
import numpy as np
import copy

def runInputOutput(param, useGUI):
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
    import FELiCS.Equation.EquationCollection as EquationCollection

    ## Initialization
    # mesh
    mesh=param.BCs.getMesh()
    # FEMSpaces
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )
    # read in mean flow
    MeanFlow = meanFlowClass(param, FEMSpaces, mesh)
    MeanFlow.importDataFromFile()
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        MeanFlow.exportBaseFlowAsHDF5()
    # export mean flow in "h5" file
    meanflowFilename = 'meanflow.h5'
    MeanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)


    ## create equations, discretize and solve
    equationColl = EquationCollection.EquationCollectionClass(
                    param,
                    FEMSpaces,
                    MeanFlow,
                    mesh
                    )
    LinearAlgebraObj = equationColl.DiscretizeFlow()

    fluctSolutList = LinearAlgebraObj.solveInputOutput()


    ## export
    if useGUI:
        ExportGUI(param, fluctSolutList, MeanFlow,FEMSpaces, equationColl,mesh)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow)
