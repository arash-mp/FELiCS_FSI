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
from os import path
import  json
from 	FELiCS.Misc.logging import Logger

# Get the logger
logger = Logger.get_logger("felics")

class MixtureClass():
    """
    Defines a thermochemical mixture for reactive flow simulations.

    This class manages properties of a mixture defined in an external JSON file,
    including species, reaction mechanisms, transport properties, and viscosity.

    **Initialize the MixtureClass object**

    Parameters
    ----------
    mixFilePath : str
        Relative path to the mixture JSON file (e.g., 'Mixture.json').
    speciesFilePath : str
        Relative path to the species JSON file (not used in the constructor but stored).

    Attributes
    ----------
    __Species__ : dict
        Dictionary of species involved in the mixture.
    __Reaction_mechanism__ : dict
        Dictionary describing the type of reaction mechanism used.
    __Pr__ : float
        Prandtl number of the mixture.
    __Viscosity__ : dict
        Dictionary containing the viscosity model and its parameters.

    Notes
    -----
    - If the mixture file is missing, default values are used.
    - A warning is issued if the mixture file is not a `.json` file.
    - The species data should be read later using `readSpeciesDict`.

    Raises
    ------
    Warning
        If the mixture file format is incorrect.
    Info
        If the file is not found and defaults are used.
    """

    def __init__(self,
    mixFilePath,
    speciesFilePath,
    ):
        """
        Initializing a mixture object based on a mixture.json file.

        Parameters
        ----------
        mixFilePath : string
            relative path to Mixture.json file
        speciesFilePath : string
            relative path to Species.json file
        """
        # set default values
        self.__Species__={}
        self.__Reaction_mechanism__ = {'type':'None'}
        self.__Pr__ = 1.0
        self.__Viscosity__ = {'type':'Constant','Constants':{'nu':1.0}}
        if path.isfile(mixFilePath):
            if mixFilePath.endswith(".json"):
                mixFile = open(
                mixFilePath,
                'r',
                )
                data = json.load(mixFile)
                for setting,value in data.items():
                    setattr(
                    self,
                    '__'+setting+'__',
                    value,
                    )
                mixFile.close()
            else:
                logger.warning('Mixture file ('+mixFilePath+') does not have the correct format and will not be read (should be a "json" file). This could lead to an unexplained error later, if the calculation is depending on data given in the mixture file. Please convert the mixture file (examples can be found in the felics-test repository).')
        else:
            logger.info('No mixture file '+mixFilePath+', using defaults.')
        
    def get_reaction_mechanism(self):
        """
        Return the dictionary describing the reaction mechanism.

        Returns
        -------
        dict
            Dictionary containing the type and details of the reaction mechanism.
        """

        return self.__Reaction_mechanism__
    
    def get_species_list(self,
    keyword='all',
    ):
        """
        Return a list of species based on a specified category.

        Parameters
        ----------
        keyword : {'all', 'transported', 'constraint'}, optional
            Filter for returned species:
            - 'all': return all species.
            - 'transported': return only species for which a transport equation is solved.
            - 'constraint': return only passive (non-transported) species.

        Returns
        -------
        list of str
            List of species names matching the specified keyword.
        """

        if keyword=='all':
            speciesList = list(self.__Species__.keys())
        elif keyword in ['transported','constraint']:
            species = self.__Species__
            speciesList=[]
            for sp in list(species.keys()):
                if species[sp]['calc']==keyword:
                    speciesList.append(sp)
        return speciesList
    
    def read_species_dict(self,
    filename,
    ):
        """
        Read species data from a dictionary file and store relevant transport properties.

        Parameters
        ----------
        filename : str
            Path to the file containing a dictionary of species and their properties.

        Notes
        -----
        The file is expected to define a dictionary with entries for each species,
        each containing keys: 'mol_weight', 'Sc', and 'Sc_t'. Only species classified as
        'transported' will be processed.
        """

        fileSpecies = open(
        filename,
        'r',
        )
        SpeciesDict = eval(fileSpecies.read())
        self.__M={}
        self.__Sc={}
        self.__Sc_t={}
        for specie in self.get_species_list('transported'):
            self.__M[specie]=SpeciesDict[specie]['mol_weight']
            self.__Sc[specie]=SpeciesDict[specie]['Sc']
            self.__Sc_t[specie]=SpeciesDict[specie]['Sc_t']
        
    
    def sc(self,
    specie,
    ):
        """
        Return the Schmidt number for a specific species.

        Parameters
        ----------
        specie : str
            Name of the species.

        Returns
        -------
        float
            Schmidt number of the species.

        Raises
        ------
        KeyError
            If the species is not defined in the species dictionary.
        """

        return self.species[specie]['Sc']

    def get_species_dict(self):
        """
        Return the dictionary of species in the mixture.

        Returns
        -------
        dict
            Dictionary with species names as keys and their properties as values.
        """
            
        return self.__Species__

    @property
    def reaction_mechanism(self):
        """
        Get the reaction mechanism definition.

        Returns
        -------
        dict
            Dictionary describing the reaction mechanism.
        """

        return self.__Reaction_mechanism__

    @property
    def pr(self):
        """
        Get the Prandtl number of the mixture.

        Returns
        -------
        float
            Prandtl number value.
        """
        
        return self.__Pr__

    @property
    def viscosity(self):
        """
        Get the viscosity model and its parameters.

        Returns
        -------
        dict
            Dictionary containing the viscosity model (e.g., constant) and its parameters.
        """
        return self.__Viscosity__

    @property
    def species(self):
        """
        Get the species dictionary (alias to `Species`).

        Returns
        -------
        dict
            Dictionary of species in the mixture.
        """

        return self.__Species__

    @property
    def species(self):
        """
        Get the species dictionary.

        Returns
        -------
        dict
            Dictionary of species in the mixture.
        """

        return self.__Species__
