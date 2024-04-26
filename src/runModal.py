import pdb
import numpy as np
import copy

def runModal(param, useGUI):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
        useGUI: Boolean, True if program is run using GUI, False if run from
        terminal directly
    '''
    import FELiCS.Import as Import
    import FELiCS.DefineFEMSpaces as DefineFEMSpaces
    from FELiCS.ExportSolution import ExportGUI,ExportFromFile
    from FELiCS.meanFlowClass import meanFlowClass
    import numpy as np
    from FELiCS.fluctuationClass import fluctuationSolutions
    import copy
    import FELiCS.WeakFormulationCollection as WeakFormulationCollection
    from FELiCS.functions import (
        printError,
        printWarning,
        printDebug,
        )

    mesh=param.BCs.getMesh()
    printDebug(True,'--------------------------------')
    printDebug(True,'-- Defining FEMSpaces...')
    FEMSpaces = DefineFEMSpaces.FEMSpacesClass(
                param,
                mesh,
                )

    printDebug(True,'--------------------------------')
    printDebug(True,'-- Reading InputFlow...')
    MeanFlow = meanFlowClass(param, FEMSpaces, mesh)
    MeanFlow.importDataFromFile()
    if not param.FlowInput.MeanFlowFilePath.split('.')[-1] == 'hdf5':
        MeanFlow.exportBaseFlowAsHDF5()

    # export the mapped Meanflow:
    meanflowFilename = 'meanflow.h5'
    MeanFlow.mapToExportMeshAndExport(FEMSpaces, meanflowFilename)

    printDebug(True,'--------------------------------')
    printDebug(True,'-- Discretizing the Equations...')
    WeakFormulation = WeakFormulationCollection.WeakFormulationCollectionClass(
                    param,
                    FEMSpaces,
                    MeanFlow,
                    mesh
                    )
    LinearAlgebraObj = WeakFormulation.DiscretizeFlow()


    fluctSolutList = LinearAlgebraObj.\
                     solveGEVP(param.Case.CalculateAdjoint)


    if useGUI:
        ExportGUI(param, fluctSolutList, MeanFlow,FEMSpaces, WeakFormulation)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow)
