from abc import ABC, abstractmethod
from h5py import File
import pdb
import json
from FELiCS.Equation.MixtureClass import MixtureClass
from FELiCS.Equation.Reactions.reactionMechanism import reactionMechanismClass

class Settings(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def getAllSettingsDict(self):
        pass

    def importSettings(self,settingFilePath):
        """
        Imports settings from .json file into SubSettingsClass

        Parameters
        ----------
        settingFilePath : txt
            path to .json file containing settings
        """
        file = open(settingFilePath)
        data = json.load(file)
        for setting,value in data.items():
            if setting == 'EigenValueGuess':    # EigenValueGuess specific handling to import complex numbers
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
            setattr(self,setting,value)
        file.close()
        self.Mixture = MixtureClass(
            self.MixtureFilePath,
            self.SpeciesFilePath,
            )