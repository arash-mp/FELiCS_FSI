from FELiCS.Misc.functions import printWarning,printOK
import json

class MixtureClass():
    '''The mixture class defines a mixture
    '''
    def __init__(self,mixFilePath,speciesFilePath):
        """
        Initializing a mixture object based on a mixture.json file.

        Parameters
        ----------
        mixFilePath : string
            relative path to Mixture.json file
        speciesFilePath : string
            relative path to Species.json file
        """
        from os import path
        # set default values
        self.__Species__={}
        self.__Reaction_mechanism__ = {'type':'None'}
        self.__Pr__ = 1.0
        self.__Viscosity__ = {'type':'Constant','Constants':{'nu':1.0}}
        if path.isfile(mixFilePath):
            mixFile = open(mixFilePath, 'r')
            data = json.load(mixFile)
            for setting,value in data.items():
                setattr(self,'__'+setting+'__',value)
            mixFile.close()
        else:
            printOK('The mixture file path ('+mixFilePath+') does not point to a mixture file!')
        
   
    def getReactionMechanism(self):
        ''' Function returning the reaction mechanism '''
        return self.__Reaction_mechanism__
    def getSpeciesList(self,keyword='all'):
        ''' Function giving a list of all species in the mixture.
        Input:
        \t keyword: either...
        \t \t ... all: All species are returned
        \t \t ... transported: Only the species for which a transport equation is solved are returned
        \t \t ... constraint: Only passive species are returned
        Output:
        \t speciesList: List of all species in the mixture which agree with the keyword
        '''
        if keyword=='all':
            speciesList = list(self.__Species__.keys())
        elif keyword in ['transported','constraint']:
            species = self.__Species__
            speciesList=[]
            for sp in list(species.keys()):
                if species[sp]['calc']==keyword:
                    speciesList.append(sp)
        return speciesList
    def readSpeciesDict(self,filename):
        fileSpecies = open(filename,'r')
        SpeciesDict = eval(fileSpecies.read())
        self.__M={}
        self.__Sc={}
        self.__Sc_t={}
        for specie in self.getSpeciesList('transported'):
            self.__M[specie]=SpeciesDict[specie]['mol_weight']
            self.__Sc[specie]=SpeciesDict[specie]['Sc']
            self.__Sc_t[specie]=SpeciesDict[specie]['Sc_t']
        
    
    def Sc(self,specie):
        return self.Species[specie]['Sc']

    def getSpeciesDict(self):
        ''' Function returning the private Species Dict
        Output:
        \t speciesDict: Dictionary with all species'''
        return self.__Species__

    @property
    def reactionMechanism(self):
        ''' Function returning the reaction mechanism '''
        return self.__Reaction_mechanism__

    @property
    def Pr(self):
        ''' Function returning the Prandtl number of the mixture '''
        return self.__Pr__

    @property
    def Viscosity(self):
        ''' Function returning the Prandtl number of the mixture '''
        return self.__Viscosity__

    @property
    def species(self):
        ''' Function returning the Species number of the mixture '''
        return self.__Species__

    @property
    def Species(self):
        ''' Function returning the Species number of the mixture '''
        return self.__Species__
