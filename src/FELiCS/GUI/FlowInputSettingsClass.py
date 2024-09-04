from FELiCS.GUI.SettingsClass import Settings
import json
class FlowInputSettingsClass(Settings):
    def __init__(self):
        '''Initializing the flow input settings class
        '''
        super().__init__()
        self._settingsKind = 'FlowInput'
        FlowInputSettingsDict=self.getAllSettingsDict()
        for key in list(FlowInputSettingsDict.keys()):
            tempStr='self.'+key+'=\"'+str(FlowInputSettingsDict[key]['default'])+'\"'
            exec(tempStr)

    def getAllSettingsDict(self):
        '''Function returning all flow input settings with default values'''
        CaseSettingsDict={\
            'MeanFlowFilePath':{'datatype':str,'default':''},\
            'AveragingDirection':{'datatype':str,'default':'None'},\
            'AveragingAxis':{'datatype':str,'default':'x'}\
        }
        return CaseSettingsDict

    def importSettings(self,settingFilePath):
        ''' Loading Mean Flow parameters from file '''
        FlowInputSettingsDict=self.getAllSettingsDict()
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
        if isfile(self.MeanFlowFilePath):
            EverythingPresent=True
        else:
            printOK('Set input flow file!')
            EverythingPresent=False
        return EverythingPresent
