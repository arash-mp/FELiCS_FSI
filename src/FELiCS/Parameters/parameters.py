#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
from h5py import File, string_dtype, AttributeManager
import pdb
import numpy as np
from FELiCS.Misc.functions import getLastGitCommit

class parameters():
    """
    Class that holds and manages the parameters for a simulation case. This includes 
    various settings related to the case, flow input, boundary conditions, IO resolvent,
    numerical methods, and export options. The class provides methods to export and 
    import settings from files, check parameter completeness, and manage older parameter formats.
    """
    def __init__(self):
        """
        Initialize the parameters object by loading the necessary settings classes 
        for case, flow input, boundary conditions, IO resolvent, numerics, and export.

        This constructor imports several classes from the GUI module and initializes 
        them as attributes of the parameters class.
        """
        from FELiCS.GUI.CaseSettingsClass import CaseSettingsClass
        self.Case=CaseSettingsClass()
        from FELiCS.GUI.FlowInputSettingsClass import FlowInputSettingsClass
        self.FlowInput=FlowInputSettingsClass()
        from FELiCS.GUI.BCsSettingsClass import BCsSettingsClass
        self.BCs = BCsSettingsClass(self.Case.MeshFilePath)
        from FELiCS.GUI.IOResolventSettingsClass import IOResolventSettingsClass
        self.IOResolvent = IOResolventSettingsClass()
        from FELiCS.GUI.NumericsSettingsClass import NumericsSettingsClass
        self.Numerics = NumericsSettingsClass()
        from FELiCS.GUI.ExportSettingsClass import ExportSettingsClass
        self.Export = ExportSettingsClass()
        pass
    def export(self, filestring):
        """
        Export the parameters to a specified file.

        Parameters
        ----------
        filestring : str
            The path of the file to export the parameters to. If the file has a 
            '.h5' extension, the parameters will be stored in an HDF5 format. 
            Otherwise, the parameters will be written as a plain text file.
        """
        from inspect import isclass

        # Writing the parameters to a file
        if '.h5' in filestring:
            file = File(filestring, 'a')

            if 'parameters' not in file.keys():
                paramGroup = file.create_group('parameters')
            else:
                paramGroup = file['parameters']
            #dt = string_dtype()
            # add the fenics-version:
            if 'FELiCSVersion' not in paramGroup.keys():
                felicsVersion = paramGroup.create_group('FELiCSVersion')
            else:
                felicsVersion = paramGroup['FELiCSVersion']
            felicsVersion.attrs.create('FELiCSVersion', data=getLastGitCommit())
        else:
            file = open(filestring,'w')

        # Iterate over all attributes of the object
        for group in dir(self):
            # Only not internal (marked by _) and not callable attributes should be exported
            if not group.startswith('_') and not callable(eval('self.'+group)):
                # Iterate over all attributes of the object
                for parameter in dir(eval('self.'+group)):
                    # Only not internal (marked by _) and not callable or class attributes should be exported
                    if not parameter.startswith('_') and not callable(eval('self.'+group+'.'+parameter)) and not isclass(parameter):

                        # Check if the attributes are string, float, int, bool or list and export them
                        if type(eval('self.'+group+'.'+parameter)) == str:
                            if isinstance(file, File):
                                if group not in list(paramGroup.keys()):
                                    currentGroup = paramGroup.create_group(group)
                                    # stringInArray = np.array([eval(f'self.{group}.{parameter}')])
                                    #

                                    currentGroup.attrs.create(parameter, data=eval(f'self.{group}.{parameter}'))
                                else:

                                    #tringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    file[f'parameters/{group}'].attrs.create(parameter, data=eval(f'self.{group}.{parameter}'))
                            else:
                                file.writelines(parameter+'='+'\''+str(eval('self.'+group+'.'+parameter))+'\''+'\n')
                        elif type(eval('self.'+group+'.'+parameter)) in [float,int,bool,list,dict]:
                            if isinstance(file, File):
                                if group not in list(paramGroup.keys()):
                                    currentGroup = paramGroup.create_group(group)
                                    #try:
                                    currentGroup.attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
                                    # except:
                                    #   # stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    #
                                    #   currentGroup.attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
                                else:
                                    #try:
                                    file[f'parameters/{group}'].attrs.create(parameter, data=str(eval(f'self.{group}.{parameter}')))
                                    # except:
                                    #   stringInArray = np.array(str(eval(f'self.{group}.{parameter}')))
                                    #   file[f'param/{group}'].create_dataset(parameter, data=stringInArray)
                            else:
                                file.writelines(parameter+'='+str(eval('self.'+group+'.'+parameter))+'\n')


        file.close()
    def importFromFile(self,filestring):
        """
        Import parameters from a specified file.

        Parameters
        ----------
        filestring : str
            The path of the parameter file to import. If the file has a '.h5' 
            extension, the parameters will be loaded from an HDF5 file. Otherwise, 
            the parameters will be imported from a plain text file.
        """
        if '.h5' in filestring:
            self.Case.importFromH5File(filestring)
            self.FlowInput.importFromH5File(filestring)
            self.BCs.importFromH5File(filestring)
            self.BCs.readDomainData(
                                self.Case.MeshFilePath,
                self.Case.getExtendedTransportedQuantityList(),
                self.Case.CoordinateSystem,
                                self.Case.m,
                                )
            self.IOResolvent.importFromH5File(filestring)
            self.Numerics.importFromH5File(filestring)
            self.Export.importFromH5File(filestring)
        else:

            self.Case.importSettings(filestring)
            self.FlowInput.importSettings(filestring)
            self.BCs.importSettings(filestring)
            self.BCs.readDomainData(
                                self.Case.MeshFilePath,
                                self.Case.nDim,
                                self.Case.getExtendedTransportedQuantityList(),
                                self.Case.CoordinateSystem,
                                self.Case.m,
                                )
            self.IOResolvent.importSettings(filestring)
            self.Numerics.importSettings(filestring)
            self.Export.importSettings(filestring)
            self.getOldParameters()
    def complete(self):
        """
        Check all parameters for completeness and consistency.

        This method checks if all the required parameters are set correctly 
        and ensures that the settings are consistent. It also checks the 
        completeness of specific sub-groups based on the analysis mode.

        Returns
        -------
        bool
            True if all parameters are complete and consistent, False otherwise.
        """
        #Check all sub parameter groups for completeness.
        CompleteBool=False
        if self.Case.complete():
            CompleteBool = self.FlowInput.complete() and\
               self.BCs.complete(self.Case.getExtendedTransportedQuantityList) and\
               self.Numerics.complete(self.Case.getTransportedQuantityList())
        # IOResolvent is only checked if Analysis mode is wither of IO or Resolvent
            if self.Case.AnalysisMode in ['Input-Output','Resolvent']:
                CompleteBool= CompleteBool and self.IOResolvent.complete()
        return CompleteBool

    def getOldParameters(self):
        """
        Provide the parameters in the 'old' format for compatibility with the rest 
        of the code. This function is deprecated and will eventually be removed.

        The function retrieves the old parameters that were set in an older fashion 
        for backward compatibility. In the future, this function will be obsolete.

        Notes
        -----
        This method is considered redundant and should be removed once the rest 
        of the code is adapted to the new parameter handling system.
        """        # AdditionalVelocityComponents
        from copy import copy
        from FELiCS.Misc.functions import printDebug, printError
        #from fenics import DirichletBC,MeshFunction
        from dolfinx.fem import dirichletbc as DirichletBC
        from FELiCS.SpaceDisc.DefineFEMSpaces import FEMSpacesClass

        #nVelocityComponents
        self.nVelocityComponents=len(self.Case.getInternalVelocityComponents())
        #SolutionList
        self.SolutionList=self.Case.getTransportedQuantityList()
        #
        self.VelocityComponents = self.Case.getInternalVelocityComponents()
        #SpeciesList
        self.SpeciesList=self.Case.Mixture.getSpeciesList('transported')


        self.ExportMode='both'

        #debug
        self.debug=True
        #NumericalScheme
        self.NumericalScheme = 'Continuous Galerkin'
        #Additional Species
        self.additionalSpecies=[]
