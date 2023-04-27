from .SettingsClass import Settings
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

			#Read whole file
			lines = file.readlines()
			# Add every line of the file as an attribute to the object
			for line in lines:
				if line.split('=')[0].strip() in list(SettingsDict.keys()):
					exec('self.'+line)
			file.close()

	def complete(self,variableList):
		''' Checking if all necessary case attributes are present '''
		from functions import printOK
		from os.path import isdir
		#Check inputs for completeness and correctness ...'
		EverythingPresent=True
		if self.Folder=='' or not isdir(self.Folder):
			EverythingPresent=False
			printOK('Set export folder!')
		return EverythingPresent
