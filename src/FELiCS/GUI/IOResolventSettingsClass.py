from FELiCS.GUI.SettingsClass import Settings
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
            'ForcingNorm':{'datatype':str,'default':'TKE'},
            'ResponseNorm':{'datatype':str,'default':'TKE'},
            'ForcingCoeff':{'datatype':list,'default':[]},
            'ResponseCoeff':{'datatype':list,'default':[]}\
        }
        return CaseSettingsDict

    def importSettings(self,settingFilePath):
        ''' Loading Mean Flow parameters from file '''
        SettingsDict=self.getAllSettingsDict()
        if not settingFilePath =='':
            file = open(settingFilePath)

            #Read whole file
            lines = file.readlines()
            # Add every line of the file as an attribute to the object
            for line in lines:
                if line.split('=')[0].strip() in list(SettingsDict.keys()):
                    exec('self.'+line)
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
            printOK('Set variables used for forcing term in ForcingCoeff')
        if not self.ResponseCoeff:
            EverythingPresent=False
            printOK('Set variables used for response term in ResponseCoeff')
        if not self.ResponseNorm:
            EverythingPresent=False
            printOK('Set norm type for response in ResponseNorm')
        if not self.ForcingNorm:
            EverythingPresent=False
            printOK('Set norm type for response in ForcingNorm')

        return EverythingPresent
