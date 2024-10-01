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
        if not settingFilePath =='':
            file = open(settingFilePath)
            data = json.load(file)
            for setting,value in data.items():
                if setting == 'EigenValueGuess' and isinstance(value, list):
                    EV_guesses=[]
                    for EV_guess in value:
                        if isinstance(EV_guess, list):
                            EV_guesses.append(complex(EV_guess[0],EV_guess[1]))
                        else:
                            EV_guesses.append(EV_guess)
                        setattr(self,setting,EV_guesses)
                else:
                    setattr(self,setting,value)
            file.close()
            #The Mixture and reaction are not loaded but constructed from the inputs
            self.Mixture = MixtureClass(
                self.MixtureFilePath,
                self.SpeciesFilePath,
                )
            #print(self.Mixture.getReactionMechanism())
            #self.reactionMechanism = reactionMechanismClass(self.Mixture.getReactionMechanism()['type'])