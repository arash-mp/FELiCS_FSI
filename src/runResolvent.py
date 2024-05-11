import pdb
import numpy as np
import copy

def runResolvent(param, useGUI):
    '''This function runs the calculations preset in param
    Input:
        param: Parameter objects (see parameters.py), defining the case
        useGUI: Boolean, True if program is run using GUI, False if run from
        terminal directly
    '''
    import FELiCS.IO.Import as Import
    import FELiCS.SpaceDisc.DefineFEMSpaces as DefineFEMSpaces
    from FELiCS.IO.ExportSolution import ExportGUI,ExportFromFile
    from FELiCS.Fields.meanFlowClass import meanFlowClass
    import numpy as np
    from FELiCS.Fields.fluctuationClass import fluctuationSolutions
    import copy
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


    fluctSolutList = LinearAlgebraObj.solveResolvent(
            equationColl)

    ## export
    if useGUI:
        ExportGUI(param, fluctSolutList, MeanFlow,FEMSpaces, equationColl,mesh)
    else:
        ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow)
