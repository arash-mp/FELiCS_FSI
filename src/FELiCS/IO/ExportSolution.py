#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * opermission of the copyright owner, the FLOW group at TU Berlin
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de
# */
#%%
'''
# **********************************************************************
# * This file deals with the Export of the Results data
# *
# * This file was created by Thomas L. Kaiser. Significant contributions
# * were made by
# * -
# *
# *
# ********************
'''

import  matplotlib
import  numpy                               as np
import  tkinter                             as tk
import  matplotlib.pyplot                   as plt
from    dolfinx                             import plot
from    dolfinx.fem                         import Function
from    matplotlib.backends.backend_tkagg   import FigureCanvasTkAgg, NavigationToolbar2Tk
from    FELiCS.Misc.functions               import getLastGitCommit
from 	FELiCS.Misc.logging                 import Logger

# Get the logger
logger = Logger.get_logger("felics")

# NOTE (Simon) is this needed?
matplotlib.use("Agg")

def writeCSVGains(param,gains):
    ''' This file writes the gains to a CSV file in the solution directory provided by the user'''
    import csv
    with open(param.Export.ExportFolder+'/gains.csv', mode='w') as writer_file:
        writer          = csv.writer(writer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        n_gains,n_omega = np.shape(gains)
        head            = ['omega']
        for i in range(n_gains):
            head.append('gain'+str(i))
        writer.writerow(head)
        for i in range(n_omega):
            line = [param.IOResolvent.Omegas[i]]
            for j in range(n_gains):
                line.append(gains[j,i])
            writer.writerow(line)
    logger.debug("Gains saved in: %s/gains.csv" % param.Export.ExportFolder)

def writeCSVSpectrum(param,spectrumDirect,spectrumAdjoint=[0]):
    ''' This file writes the gains to a CSV file in the solution directory provided by the user'''
    import csv
    #TODO Sophie: how to handle different numbers of solutions for adjoint/direct?
    if param.Case.CalculateAdjoint and len(spectrumDirect)==len(spectrumAdjoint):

        with open(param.Export.ExportFolder+'/spectrum.csv', mode='w') as writer_file:
            writer = csv.writer(writer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            head=['omega_direct_r','omega_direct_i','omega_adjoint_r','omega_adjonint_i']
            writer.writerow(head)
            for i in range(0,len(spectrumDirect),len(spectrumAdjoint)):
                line=[str(np.real(spectrumDirect[i])),str(np.imag(spectrumDirect[i])),str(np.real(spectrumAdjoint[i])),str(np.imag(spectrumAdjoint[i]))]
                writer.writerow(line)

    else:
        with open(param.Export.ExportFolder+'/spectrum.csv', mode='w') as writer_file:
            writer = csv.writer(writer_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            head=['omega_direct_r','omega_direct_i']
            writer.writerow(head)
            for i in range(0,len(spectrumDirect)):
                line=[str(np.real(spectrumDirect[i])),str(np.imag(spectrumDirect[i]))]
                writer.writerow(line)
    logger.debug("Spectrum saved in: %s/spectrum.csv" % param.Export.ExportFolder)


def writeLastGitCommit(param):
    commitLabel=getLastGitCommit()
    gitFile=open(param.Export.ExportFolder+'/LastCommit.txt', mode='w')
    gitFile.write(commitLabel)


def ExportFromFile(param,FEMSpaces,fluctSolutList,MeanFlow):
    # automated export from here on-Felics call wihtout gui --> no export-window: all frequencies and the first 2 leading modes
    # xy = FEMSpaces.exportMesh.coordinates()
    # mesh_cells= FEMSpaces.exportMesh.cells()
    # in the case of resolvent analysis: expport all frequencies and the first two leading modes

    if param.Case.AnalysisMode in ['Resolvent','Input-Output']:
        # get the gains from the list of fluctuation Solutions:
        # gains   = np.zeros((param.Numerics.nSolut, len(param.IOResolvent.Omegas)))
        # NOTE: Removed the gains export to file from here, as it is now done by modeCollection.exportSpectrumToCSV
        ct      = 1
        for fluctSolut in fluctSolutList:
            logger.debug(f"Saving output file {ct:0.0f} / {len(fluctSolutList):0.0f}")
            # if fluctSolut.solutionKind == 'Response':
            #     indexOmega = param.IOResolvent.Omegas.index(fluctSolut.omega)
            #     gains[fluctSolut.gainNumber, indexOmega] = fluctSolut.gainValue
            fileName = f'{param.Case.AnalysisMode}_Omega{np.round(fluctSolut.omega, 3)}_{fluctSolut.solutionKind}_gain{fluctSolut.gainNumber}.h5'
            fluctSolut.exportSolution(fileName, 'o')
            ct += 1
        # writeCSVGains(param,gains)
        writeLastGitCommit(param)

    # in case of modal analysis: export the whole spectrum and every mode
    elif param.Case.AnalysisMode == 'Modal':

        # NOTE: Removed the gains export to file from here, as it is now done by modeCollection.exportSpectrumToCSV
        # eValDirect = []
        # eValAdjoint = []

        # for fluctSolut in fluctSolutList:
        #     if fluctSolut.solutionKind == 'Direct':
        #         eValDirect.append(fluctSolut.omega)
        #     else:
        #         eValAdjoint.append(fluctSolut.omega)

        # convert the list of eVals to a np-array:
        # eValDirect      = np.array(eValDirect)
        # idxDirect       = np.argmax(np.imag(eValDirect))
        
        # if param.Case.CalculateAdjoint:
        #     eValAdjoint = np.array(eValAdjoint)
        #     # idxAdjoint  = np.argmax(np.imag(eValAdjoint))
        #     writeCSVSpectrum(param,eValDirect,eValAdjoint)
        # else:
        #     writeCSVSpectrum(param,eValDirect)  

        for fluctSolut in fluctSolutList:
            fileName    = f'{param.Case.AnalysisMode}Solution_Omega_{fluctSolut.solutionKind}_{np.round(fluctSolut.omega, 3)}.h5'
            fluctSolut.exportSolution(fileName, 'o')
        writeLastGitCommit(param)
