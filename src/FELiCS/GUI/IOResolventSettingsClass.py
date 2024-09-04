from FELiCS.GUI.SettingsClass import Settings
import json
class IOResolventSettingsClass(Settings):
    def __init__(self):
        '''Initializing the input-output/resolvent settings class
        '''
        super().__init__()
        self._settingsKind = 'IOResolvent'      
        SettingsDict=self.getAllSettingsDict()
        for key in list(SettingsDict.keys()):
            if type(SettingsDict[key]['default']) ==list:
                quotation=''
            else:
                quotation='\"'
            tempStr='self.'+key+'='+quotation+str(SettingsDict[key]['default'])+quotation
            exec(tempStr)

    def getAllSettingsDict(self):
        '''Function returning all input-output/resolvent settings with default values'''
        CaseSettingsDict={\
            'Omegas':{'datatype':list,'default':[]},
            'ForcingBoundaryIndices':{'datatype':list,'default':[]},
            'ForcingMode':{'datatype':str,'default':'Body'},
            'ForcingCoeff':{'datatype':list,'default':[]},
            'ResponseCoeff':{'datatype':list,'default':[]}\
        }
        return CaseSettingsDict

    def importSettings(self,settingFilePath):
        ''' Loading Mean Flow parameters from file '''
        SettingsDict=self.getAllSettingsDict()
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

    def complete(self):
        ''' Checking if all necessary case attributes are present '''
        from os.path import isfile
        from FELiCS.Misc.functions import printOK
        #Only the mesh is absolutely necessary...'
        EverythingPresent=True
        if not len(self.Omegas)>0:
            EverythingPresent=False
            printOK('Set frequencies!')
        if not self.ForcingCoeff:
            EverythingPresent=False
            printOK('Set forcing norms!')
        if not self.ResponseCoeff:
            EverythingPresent=False
            printOK('Set response norms!')

        return EverythingPresent
