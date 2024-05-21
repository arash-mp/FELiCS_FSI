from FELiCS.GUI.SettingsClass import Settings
from FELiCS.Equation.MixtureClass import MixtureClass
from FELiCS.Equation.Reactions.reactionMechanism import reactionMechanismClass
class CaseSettingsClass(Settings):
    def __init__(self):
        '''Initializing the Case settings '''
        # First initialize parameters to default
        # _CaseSettingsDict_ contains all parameters and references their default value
        super().__init__()
        self._settingsKind = 'Case'
        CaseSettingsDict=self.getAllSettingsDict()
        for key in list(CaseSettingsDict.keys()):
            if type(CaseSettingsDict[key]['default']) == str:
                tempStr='self.'+key+'=\"'+str(CaseSettingsDict[key]['default'])+'\"'
            else:
                tempStr='self.'+key+'='+str(CaseSettingsDict[key]['default'])
            exec(tempStr)
        self.reactionMechanism = reactionMechanismClass('None') 
        self.customSolutions=[]


    def getAllSettingsDict(self):
        '''Function returning all CaseSettings with default values'''
        CaseSettingsDict={\
            'nDim':{'datatype':int,'default':2},
            'CoordinateSystem':{'datatype':str,'default':'Cartesian'},
            'TurbulenceModel':{'datatype':str,'default':'None'},
            'm':{'datatype':int,'default':0},
            'MolViscModel':{'datatype':str,'default':'Constant'},
            'MolViscPerturbModel':{'datatype':str,'default':'None'},
            'MolVisc':{'datatype':int,'default':0.0},
            'AnalysisMode':{'datatype':str,'default':'Modal'},
            'CalculateAdjoint':{'datatype':bool,'default':True},
            'TransVelFluc':{'datatype':bool,'default':False},
            'SetOfEquations':{'datatype':dict,'default':{
                'Momentum':{'Equation':'NSPrimitive','Variable':'u'},
                'Mass': {'Equation':'Continuity','Variable':'p'},
                'Energy'    : {'Equation':'None','Variable':'None'},
                'Species': {'Equation':'None','Variable':'None'},
                'EquationOfState': {'Equation':'None'    ,'Variable':'None'},
                'Custom1': {'Equation':'None'    ,'Variable':'None'},
                'Custom2': {'Equation':'None'    ,'Variable':'None'}
                }
                },
            'Reaction':{'datatype':bool,'default':False},
            'MixtureFilePath':{'datatype':str,'default':''},
            'SpeciesFilePath':{'datatype':str,'default':''},
            'MeshFilePath':{'datatype':str,'default':''}
        }
        return CaseSettingsDict

            #'SetOfEquations':{'datatype':dict,'default':{'Navier-Stokes':'Primitive Variables','Energy':'None','Species': 'None', 'equationOfState': 'None'}}
    def importSettings(self,settingFilePath):
        ''' Loading Case parameters from file '''
        from FELiCS.Equation.Reactions.reactionMechanism import reactionMechanismClass
        CaseSettingsDict=self.getAllSettingsDict()
        if not settingFilePath =='':
            file = open(settingFilePath)

            #Read whole file
            lines = file.readlines()
            #Add every line of the file as an attribute to the object
            for line in lines:
                if line.split('=')[0].strip() in list(CaseSettingsDict.keys()):
                    exec('self.'+line)
            file.close()
            #The Mixture and reaction are not loaded but constructed from the inputs
            self.Mixture = MixtureClass(
                self.MixtureFilePath,
                self.SpeciesFilePath,
                )
            #print(self.Mixture.getReactionMechanism())
            #self.reactionMechanism = reactionMechanismClass(self.Mixture.getReactionMechanism()['type']) 
    def importFromH5File(self, h5FileName):
        """
        This function adapts the inherited function of same name from the SettingsClass
        """
        from FELiCS.Equation.Reactions.reactionMechanism import reactionMechanismClass
        super().importFromH5File(h5FileName);
        #The Mixture is not loaded but constructed from the inputs
        self.Mixture = MixtureClass(
            self.MixtureFilePath,
            self.SpeciesFilePath,
            )
        self.reactionMechanism = reactionMechanismClass(self.Mixture.getReactionMechanism()['type']) 

    def getInternalVelocityComponents(self):
        ''' Provides a list of velocity components, which are directed within the dimensions of the mesh '''
        if self.CoordinateSystem=='Cartesian':
            VelCompList = ['x','y']
            if self.nDim>2:
                VelCompList.append('z')
        elif self.CoordinateSystem=='Cylindrical':
            VelCompList = ['x','r']

        return VelCompList

    def getExternalVelocityComponents(self):
        ''' Provides a list of velocity components, which are directed outside the dimensions of the mesh '''
        if self.CoordinateSystem=='Cartesian':
            if self.m!=0 and self.nDim==2:
                VelCompList = ['z']
            else:
                VelCompList = []
        elif self.CoordinateSystem=='Cylindrical':
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
        if not self.SetOfEquations['Momentum']['Variable'] == 'None':
            SolutionList.append(self.SetOfEquations['Momentum']['Variable'])
        if not self.SetOfEquations['Mass']['Variable'] == 'None':
            SolutionList.append(self.SetOfEquations['Mass']['Variable'])
        if not self.SetOfEquations['Energy']['Variable'] == 'None':
            SolutionList.append(self.SetOfEquations['Energy']['Variable'])
        if 'Species' in list(self.SetOfEquations.keys()):
            if self.SetOfEquations['Species']['Variable'] == 'Y':
                for species in list(self.Mixture.getSpeciesList('transported')):
                    SolutionList.append(species)

        # Add custom variables if they are given in the parameter file. 
        # TODO Sophie: think of something better, also: variable number of additional equations
        try:
            if not self.SetOfEquations['Custom1']['Variable'] == 'None':
                SolutionList.append(self.SetOfEquations['Custom1']['Variable'])
        except:
            pass
        try:
            if not self.SetOfEquations['Custom2']['Variable'] == 'None':
                SolutionList.append(self.SetOfEquations['Custom2']['Variable'])
        except:
            pass
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
        ''' This function provides the mean fields which mus be read in.'''

        MeanList=[]
        # Add velocity components
        MeanList.append('u')

        # If necessary, add density and enthalpy diffusion
        if 'rho' in self.getTransportedQuantityList():
            MeanList.append('rho')
            if self.MolViscModel == 'File' or self.MolViscPerturbModel == 'Sutherland mean':
                MeanList.append('alpha')
        # Add species which are transported
        for specie in self.Mixture.getSpeciesList('transported'):
            MeanList.append(specie)
            if self.MolViscModel == 'File' or self.MolViscPerturbModel == 'Sutherland mean':
                MeanList.append('D_'+specie)

        # If Input-Output analysis is used, the forcing must be read in (at least curently) for
        # every conservative variable ()...
        if self.AnalysisMode in ['Input-Output']:
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
        if self.AnalysisMode in ['Resolvent']:
            MeanList.append('responseDomain')
            MeanList.append('forcingDomain')

        # Always look for a sponge variable in the mean flow file
        # if not present it will be zero
        MeanList.append('spg')

        # Add species, which are not transported
        for specie in self.Mixture.getSpeciesList('constraint'):
            MeanList.append(specie)

        if self.TurbulenceModel in ['File']:
            MeanList.append('nuturb')
        if self.MolViscModel in ['File'] or self.MolViscPerturbModel in ['Sutherland mean']:
            MeanList.append('nulam')
        if self.TurbulenceModel in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cylindrical':
            MeanList.extend(['rstxx', 'rstrr', 'rsttt', 'rstxr', 'rstxt', 'rstrt','rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
        elif self.TurbulenceModel in ['Boussinesq', 'TKE-based', 'Boussinesq(xr)'] and self.CoordinateSystem == 'Cartesian':
            MeanList.extend(['rstxx', 'rstyy', 'rstzz', 'rstxy', 'rstxz', 'rstyz'])
    #    if self.SetOfEquations['Energy']['Equation'] == 'Enthalpy':
    #        MeanList.append('cp')
    #        MeanList.append('alpha')
    #        MeanList.append('he')
    #        MeanList.append('T')
    #        MeanList.append('molarMass')
        if self.Reaction:
            MeanList.append('dQ')
        return MeanList

    def getNVelocityComponents(self):
        return len(self.getVelocityComponents())

    def complete(self):
        ''' Checking if all necessary case attributes are present '''
        from os.path import isfile
        from FELiCS.Misc.functions import printOK
        #Only the mesh is absolutely necessary...'
        EverythingPresent=True
        if not isfile(self.MeshFilePath):
            printOK('Set mesh file!')
            EverythingPresent=False
        if not isfile(self.MixtureFilePath):
            printOK('Set mixture file!')
        return EverythingPresent
