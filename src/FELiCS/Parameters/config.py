from os import system, path
from abc import ABC, abstractmethod
from h5py import File
import pdb
import json
from FELiCS.Equation.MixtureClass import MixtureClass
from FELiCS.Equation.Reactions.reactionMechanism import reactionMechanismClass
from FELiCS.SpaceDisc.FELiCSMesh import FELiCSMesh
from FELiCS.Misc.functions import(
    getLastGitCommit,
    printWarning
    )

class dotdict(dict):
    """
    Adds possibility to use dot notation to access entries of dictionary dict
    also works for nested dictionaries

    Parameters
    ----------
    dict : dictionary
        the dictionary which is going to become a dot-dictionary
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

class config(ABC):
    def __init__(self):
        self.__BCIDs__=[]   # move to BC
        SettingsDict = self.getAllSettingsDict()

        # Add all default fields & subfields to self
        for category in SettingsDict:
            dict = dotdict()
            for parameter in SettingsDict[category]:
                dict[parameter]=SettingsDict[category][parameter]["default"]
            setattr(self,category,dict)


    def importFromFile(self,configFilePath):
        """
        Imports parameters from .json file

        Parameters
        ----------
        configFilePath : txt
            path to .json file containing parameters
        """

        # TODO: check if we have a .json or .h5
        # LOAD data accordingly
        file = open(configFilePath)
        data = json.load(file)

        # TODO: put the rest of this in a validate_config method
        #   1: Check if we have a modal or resolvent/input-output type of analysis
        #   2: Load either "eigenvalueguess" or "omega" depending on the above
        #   3: Call a method that checks if files and dir needed do exist
        #       def check_exist_files_dirs(list_files, list_dirs):
        #   4: Loop over the default config (self) and check if fields are in data
        #       if missing, depending on field:
        #           i. set the default and print warning
        #           ii. error-print + shut down
        #   5: call a method to set calculated settings (previously from case)

        # Loop over the default dictionnary and overwrite with values from data
        for category,parameters in data.items():
            subcategory = eval("self."+category)
            for parameter,value in parameters.items():
                if parameter == 'EigenValueGuess':    # EigenValueGuess specific handling to import complex numbers
                    if isinstance(value, list):
                        if not value:
                            value = []
                        elif isinstance(value[0],str):
                            value = list(map(complex,value))
                        else:
                            value = value
                    else:
                        if isinstance(value,str):
                            value = [complex(value)]
                        else:
                            value = [value]
                subcategory[parameter] = value

                # TODO: the following in a dedicated method
                if parameter == 'BCsFilePath':
                    if not path.isfile(value):
                        printWarning("BC file not found!")
                if parameter == 'MeshFilePath':
                    if not path.isfile(value):
                        printWarning("mesh file not found!")
                if parameter == 'MeanFlowFilePath':
                    if not path.isfile(value):
                        printWarning("meanFlow file not found!")
                if parameter == 'ExportFolder':
                    if not path.isdir(value):
                        printWarning("export folder not found!")
            subcategory = dotdict(subcategory)
        file.close()

        self.Mixture = MixtureClass(
            self.Case["MixtureFilePath"],
            self.Case["SpeciesFilePath"],
            )
        
        self.readDomainData(
            self.Case["MeshFilePath"],
            self.Case["nDim"],
            self.getExtendedTransportedQuantityList(),
            self.Case["CoordinateSystem"],
            self.Case["m"]
            )
        
        self.debug = True # specify here if printDebug messages should be shown
        
        # calculated parameters: TODO: put in a method
        self.BoundaryCondition.nVelocityComponents  = len(self.getVelocityComponents())
        self.Case.SolutionList                      = self.getTransportedQuantityList()
        self.BoundaryCondition.VelocityComponents   = self.getVelocityComponents()
        # self.Mixture.SpeciesList                  = self.Mixture.getSpeciesList('transported') # not used, move to mixture
        # self.reactionMechanism                    = reactionMechanismClass(self.Mixture.getReactionMechanism()['type']) 

        # TODO: set in defaults
        self.Numerics.NumericalScheme               = 'Continuous Galerkin' 
        

    def importFromH5File(self, h5FileName):
        # TODO: load the parameters into a "data" dictionnary
        # similar to what we get from loading a .json

        hf = File(h5FileName, 'r')
        if 'param' not in hf.keys():
            return None
        settingsDict = self.getAllSettingsDict()
        if self._settingsKind in hf['param'].keys():
            for settingsParameter in list(settingsDict.keys()):
                if settingsParameter in hf[f'param/{self._settingsKind}'].attrs.keys():
                    groupName = f'param/{self._settingsKind}'
                    datatype = settingsDict[settingsParameter]['datatype']
                    value = hf[groupName].attrs[settingsParameter]
                    #pdb.set_trace()
                    if 'int' in str(datatype):
                        try:
                            value = int(float(value))
                            exec(f'self.{settingsParameter} = {value}')
                        except:
                            pdb.set_trace()
                    elif 'str' in str(datatype):
                        exec(f'self.{settingsParameter} = "{value}"')
        hf.close()

        # This will not be needed anymore
        self.Mixture = MixtureClass(
            self.Case.mixtureFilePath,
            self.SpeciesFilePath,
            )


    def export(self, filestring):
        ## __________currently not working__________
        '''function exporting the parameters to a file
        \t Input:
        \t -filestring: path of parameter file'''
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


    # @abstractmethod
    def getAllSettingsDict(self):
        '''Function returning all boundary condition settings with default values'''
        SettingsDict={
            'BoundaryCondition':{
                'BCsFilePath':          {'datatype':str,    'default':''},
                'nVelocityComponents':  {'datatype':int,    'default':2},
                'VelocityComponents':   {'datatype':list,   'default':[]}
            },
            'Case':{
                'AnalysisMode':         {'datatype':str,    'default':'Modal'},
                'CalculateAdjoint':     {'datatype':bool,   'default':True},
                'CoordinateSystem':     {'datatype':str,    'default':'Cartesian'},
                'm':                    {'datatype':int,    'default':0},
                'MeshFilePath':         {'datatype':str,    'default':''},
                'MixtureFilePath':      {'datatype':str,    'default':''},
                'MolVisc':{'datatype':int,'default':0.0},
                'MolViscModel':{'datatype':str,'default':'Constant'},
                'nDim':{'datatype':int,'default':2},
                'Reaction':{'datatype':bool,'default':False},
                'SetOfEquations':{
                    'datatype':dict,
                    'default':{
                        'Momentum': {'Equation':'NSPrimitive','Variable':'u'},
                        'Mass':     {'Equation':'Continuity','Variable':'p'},
                        'Energy':   {'Equation':'None','Variable':'None'},
                        'Species':  {'Equation':'None','Variable':'None'},
                        'EquationOfState': {'Equation':'None','Variable':'None'}
                    }
                },
                'SolutionList':{'datatype':list,'default':[]},
                'SpeciesFilePath':{'datatype':str,'default':''},
                'TransVelFluc':{'datatype':bool,'default':False},
                'TurbulenceModel':{'datatype':str,'default':'None'}
            },
            'Export':{
                'ExportFolder':{'datatype':str,'default':''},
                'hdf':{'datatype':bool,'default':False},
                'mat':{'datatype':bool,'default':False},
                'vtk':{'datatype':bool,'default':True}
            },
            'FlowInput':{
                'AveragingAxis':{'datatype':str,'default':'x'},
                'AveragingDirection':{'datatype':str,'default':'None'},
                'MeanFlowFilePath':{'datatype':str,'default':''},
            },
            'IOResolvent':{
                'ForcingBoundaryIndices':{'datatype':list,'default':[]},
                'ForcingCoeff':{'datatype':list,'default':[]},
                'ForcingMode':{'datatype':str,'default':'Body'},
                'Omegas':{'datatype':list,'default':[]},
                'ResponseCoeff':{'datatype':list,'default':[]}
            },
            'Numerics':{
                'EigenValueGuess':{'datatype':list,'default':[1.0]},
                'nCPU':{'datatype':int,'default':1},
                'nSolut':{'datatype':int,'default':3},
                'PolynomialOrder':{'datatype':dict,'default':{'u':'2'},'options':[1,2]},
                'NumericalScheme': {'datatype':str,'default':'Continuous Galerkin'}
            }
        # 'SpeciesList':{'datatype:list,default:[]}
        # 'h5':{'datatype':bool,'default':True}
        # 'dim':{'datatype':int,'default':0}
        # 'MolViscPerturbModel':{'datatype':str,'default':'None'}
        # 'Video':{'datatype':bool,'default':False}
        # 'MeanFlowFilePath':{'datatype':str,'default':''}
        # 'AveragingDirection':{'datatype':str,'default':'None'}
        # 'AveragingAxis':{'datatype':str,'default':'x'}
        # 'LinearAlgebraSolver':{'datatype':str,'default':'SLEPc'}
        # 'Preconditioner':{'datatype':str,'default':'None'}
        # 'Velfluc':{'datatype:bool','default':True}
        }
        return SettingsDict

    def initBCsDict(self,VariableList):
        ''' Initialize BCsDict '''
        from FELiCS.Misc.functions import printWarning
        BCIDList=self.__BCIDs__
        #First define local BCsDict
        BCsDict={}
        for Variable in VariableList:
            BCsDict[Variable]=[]
            for BCID in BCIDList:
                BCsDict[Variable].append({'ID':BCID,'type':'Neumann','value':0.0})
        self.__BCsDict__=BCsDict

    def importBCsDict(self,VariableList):
        ''' Import a boundary condition file with checking the consistency of BCs and mesh. To read the BCs without checking use importSettings()'''
        from FELiCS.Misc.functions import printWarning
        BCIDList=self.__BCIDs__
        #First define local BCsDict
        BCsDict={}
        filepath=self.BoundaryCondition.BCsFilePath
        if not filepath == '' and path.isfile(filepath):
            self.BoundaryCondition.BCsFilePath = filepath
            for Variable in VariableList:
                BCsDict[Variable]=[]
                for BCID in BCIDList:
                    BCsDict[Variable].append({'ID':BCID,'type':'Neumann','value':0.0})
            # Read BCFile
            BCFile=open(self.BoundaryCondition.BCsFilePath)
            importDict=json.load(BCFile)
            # Loop over all variables and IDs and if needed values present in BCFile, copy the contents to the local BCsDict
            for Variable in VariableList:
                if Variable in list(importDict.keys()):
                    BCsDict[Variable]=[]
                    for BC in importDict[Variable]:
                        if BC['ID'] in BCIDList:
                            BCsDict[Variable].append(BC)
                        else:
                            printWarning('Boundary condition of variable '+Variable+' for boundary with ID '+str(BC['ID'])+' not found in file. Choosing homogeneous Neumann instead.')

                else:
                    printWarning('Boundary conditions for variable '+Variable+' not found in file. Choosing homogeneous Neumann instead.')
            # Finally, copy local BCsDict to the object
            self.__BCsDict__=BCsDict


    def setBC(self,field,BoundaryID,BCType,BCvalue):
        ''' Setting the Boundary condition of a single variable '''
        for BC in self.__BCsDict__[field]:
            if BC['ID']== BoundaryID:
                self.__BCsDict__[field][BoundaryID]['type']=BCType
                self.__BCsDict__[field][BoundaryID]['value']=BCvalue

    
    def readBCInfo(self,MeshFilePath, felicsMesh):
        ''' Input: - MeshFilePath
        This function reads both the IDs of the boundary conditions from the mesh and stores them
        in a private list of the class and also the boundary nodes and stores them in __boundaries__'''
        from numpy import unique
        # get a list of all kinds of BC indices
        self.__BCIDs__ = unique(felicsMesh.facet_tags.values)
        self.__boundaries__ = felicsMesh.facet_tags

    def complete(self,BCVariableList):
        # add completeness check later!
        return EverythingPresent

    def readDomainData(self,Meshfile,gDim,ExtendedTransportedQuantityList,coordinateSystem,m):
        ''' Input: Mesfile
        Read all the domain data from the meshfile '''
        self.readMesh(Meshfile,gDim,coordinateSystem,m)
        self.readBCInfo(Meshfile, self.__mesh__)
        self.initBCsDict(ExtendedTransportedQuantityList)
        self.importBCsDict(ExtendedTransportedQuantityList)

    def getMesh(self):
        ''' Function is returning the mesh '''
        return self.__mesh__

    def readMesh(self,MeshFile,dim,coordinateSystem,m):
        '''
        Reading Meshfile and saving it as private object

        Function Arguments:
        - MeshFile: File of a gmsh-meshfile. File needs to be in .msh format
        - gdim: Geometrical Dimension of the mesh. This argument is needed
        by the gmsh helper-functions, which read in the mesh

        Function returns:
        '''
        if not MeshFile == '' and path.isfile(MeshFile):
            self.__mesh__ = FELiCSMesh(coordinateSystem,MeshFile,dim,m)
            self.dim = self.__mesh__.gdim
        
    def getInternalVelocityComponents(self):
        ''' Provides a list of velocity components, which are directed within the dimensions of the mesh '''
        if self.Case["CoordinateSystem"]=='Cartesian':
            VelCompList = ['x','y']
            if self.Case.nDim>2:
                VelCompList.append('z')
        elif self.Case["CoordinateSystem"]=='Cylindrical':
            VelCompList = ['x','r']
        return VelCompList

    def getExternalVelocityComponents(self):
        ''' Provides a list of velocity components, which are directed outside the dimensions of the mesh '''
        if self.Case["CoordinateSystem"]=='Cartesian':
            if self.Case.m!=0 and self.Case.nDim==2:
                VelCompList = ['z']
            else:
                VelCompList = []
        elif self.Case["CoordinateSystem"]=='Cylindrical':
            VelCompList = ['t']

        return VelCompList

    def getVelocityComponents(self):
        ''' Provides a list of all velocity components, both mesh internal and external '''
        templist=self.getInternalVelocityComponents()
        templist.extend(self.getExternalVelocityComponents())

        return templist

    def getTransportedQuantityList(self):
        ''' Provides a list of all transported quantities for the given case settings '''
        SolutionList = []
        if not self.Case["SetOfEquations"]['Momentum']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Momentum']['Variable'])
        if not self.Case["SetOfEquations"]['Mass']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Mass']['Variable'])
        if not self.Case["SetOfEquations"]['Energy']['Variable'] == 'None':
            SolutionList.append(self.Case["SetOfEquations"]['Energy']['Variable'])
        if 'Species' in list(self.Case["SetOfEquations"].keys()):
            if self.Case["SetOfEquations"]['Species']['Variable'] == 'Y':
                for species in list(self.Mixture.getSpeciesList('transported')):
                    SolutionList.append(species)
        return SolutionList

    def getExtendedTransportedQuantityList(self):
        ''' Like getTransportedQuantitiyList but with all velocity components '''
        SolutionList = []
        # First add all velocity components
        transportedQuantities = self.getTransportedQuantityList()
        if 'u' in transportedQuantities:
            for component in self.getVelocityComponents():
                SolutionList.append('u'+component)
        #Then extend the list by the transported quantity list
        SolutionList.extend(self.getTransportedQuantityList())
        # Finally, remove the component 'u' if present
        if 'u' in SolutionList:
            SolutionList.remove('u')
        return SolutionList

    def getMeanFlowFieldNames(self):
        from FELiCS.Misc.functions import printDebug
        ''' This function provides the mean fields which must be read in.'''
        MeanList=[]
        # Add velocity components
        MeanList.append('u')
        # If necessary, add density and enthalpy diffusion
        if 'rho' in self.getTransportedQuantityList():
            MeanList.append('rho')
            if self.Case.MolViscModel == 'File' or self.Case.molViscPerturbModel == 'Sutherland mean':
                MeanList.append('alpha')
        # Add species which are transported
        for specie in self.Mixture.getSpeciesList('transported'):
            MeanList.append(specie)
            if self.Case.MolViscModel == 'File' or self.Case.molViscPerturbModel == 'Sutherland mean':
                MeanList.append('D_'+specie)
        # If Input-Output analysis is used, the forcing must be read in (at least curently) for
        # every conservative variable ()...
        if self.Case.AnalysisMode in ['Input-Output']:
            CurrentList=MeanList.copy()
            for entry in CurrentList:
                # If the field is the velocity vector, all components must be considered
                if entry == 'u':
                    #for component in self.VelocityComponents:
                    MeanList.append(entry+'_forcing_r')#jvs
                    MeanList.append(entry+'_forcing_i')
                else:
                    MeanList.append(entry+'_forcing_r')
                    MeanList.append(entry+'_forcing_i')
        # If Resolvent analysis is used, we always read forcing and response domains
        if self.Case.AnalysisMode in ['Resolvent']:
            MeanList.append('responseDomain')
            MeanList.append('forcingDomain')
        # Always look for a sponge variable in the mean flow file
        # if not present it will be zero
        MeanList.append('spg')
        # Add species, which are not transported
        for specie in self.Mixture.getSpeciesList('constraint'):
            MeanList.append(specie)
        if self.Case.TurbulenceModel in ['File']:
            MeanList.append('nuturb')
        if self.Case.MolViscModel in ['File'] or self.Case.molViscPerturbModel in ['Sutherland mean']:
            MeanList.append('nulam')
        if self.Case.TurbulenceModel in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.Case["CoordinateSystem"] == 'Cylindrical':
            MeanList.extend(['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt', 'rstrt','rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
        elif self.Case.TurbulenceModel in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.Case["CoordinateSystem"] == 'Cartesian':
            MeanList.extend(['rstxx', 'rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
        if self.Case.Reaction:
            MeanList.append('dQ')
        return MeanList

    def getNVelocityComponents(self):
        return len(self.getVelocityComponents())