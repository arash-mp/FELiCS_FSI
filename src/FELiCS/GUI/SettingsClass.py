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
        if not settingFilePath =='':
            file = open(settingFilePath)
            data = json.load(file)
            for key,item in data.items():
                if key == 'EigenValueGuess' and isinstance(item, list):
                    for EV_guess_raw in item:
                        if isinstance(EV_guess_raw, list):
                            EV_guess = [complex(EV_guess_raw[0],EV_guess_raw[1])]
                        else:
                            EV_guess = [EV_guess_raw]
                        setattr(self,key,EV_guess) # set attribute as list of all EV_guesses later
                else:
                    setattr(self,key,item)
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
        """

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
