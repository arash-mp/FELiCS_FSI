from FELiCS.GUI.SettingsClass import Settings
import json
class ExportSettingsClass(Settings):
    def __init__(self):
        '''Initializing the export settings class
        '''
        super().__init__()
        self._settingsKind = 'Export'
        SettingsDict=self.getAllSettingsDict()
        for key in list(SettingsDict.keys()):

            if type(SettingsDict[key]['default']) in [str]:
                quotation='\"'
            else:
                quotation=''
            tempStr='self.'+key+'='+quotation+str(SettingsDict[key]['default'])+quotation
            exec(tempStr)

    def getAllSettingsDict(self):
        '''Function returning all export settings with default values'''
        CaseSettingsDict={\
            'ExportFolder':{'datatype':str,'default':''},
            'Video':{'datatype':bool,'default':False},
            'vtk':{'datatype':bool,'default':True},
            'mat':{'datatype':bool,'default':False},
            'hdf':{'datatype':bool,'default':False},
            'h5':{'datatype':bool,'default':True}
        }
        return CaseSettingsDict

    def importSettings(self,settingFilePath):
        ''' Loading Export parameters from file '''
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

    def complete(self,variableList):
        ''' Checking if all necessary case attributes are present '''
        from FELiCS.Misc.functions import printOK
        from os.path import isdir
        #Check inputs for completeness and correctness ...'
        EverythingPresent=True
        if self.Folder=='' or not isdir(self.Folder):
            EverythingPresent=False
            printOK('Set export folder!')
        return EverythingPresent
