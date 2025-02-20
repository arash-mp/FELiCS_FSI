from FELiCS.Misc.functions import printWarning,printOK
class MixtureClass():
    '''The mixtre class defines a mixture
    '''
    def __init__(self,mixFilePath,speciesFilePath):
        '''Initializing the mixture function based on a mixture (*.mix) file.
        \t Input: mixFilePath: Providing the file path of the mixture'''
        from os import path
        #Set default values
        self.__Species__={}
        self.__Reaction_mechanism__ = {'type':'None'}
        self.__Pr__ = 1.0
        self.__Viscosity__ = {'type':'Constant','Constants':{'nu':1.0}}
        if mixFilePath == '':
            print('No mixture file chosen!'  )
        elif not path.isfile(mixFilePath):
            printOK('The mixture file path ('+mixFilePath+') does not point to a mixture file!')
        else:
            mixFile = open(mixFilePath, 'r')
            keywords=['Species','Reaction_mechanism','Pr','Viscosity']
            lineNumber=0
            while True:
                lineNumber+=1
                line = mixFile.readline()
                if not line:
                    break
                try: 
                    keyword=line.split('=')[0].strip()
                    if keyword in keywords:
                        pos=len(keyword)
                        tempstring='self.__'+line[:pos]+'__'+line[pos:]
                        exec(tempstring)    
                    else:
                        printWarning('Cannot read line '+ str(lineNumber) + ' of Mixture File. The line is ignore    d...')
                except:
                    printWarning('Cannot read line '+ str(lineNumber) + ' of Mixture File. The line is ignored...')
        #if self.getSpeciesList('transported'):
        #    self.readSpeciesDict(speciesFilePath)
   


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
