from abc import ABC, abstractmethod
from h5py import File
import pdb

class Settings(ABC):
    def __init__(self):
        pass
    @abstractmethod
    def getAllSettingsDict(self):
        pass

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
