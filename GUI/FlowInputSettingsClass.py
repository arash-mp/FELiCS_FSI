from .SettingsClass import Settings
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

			#Read whole file
			lines = file.readlines()
			# Add every line of the file as an attribute to the object
			for line in lines:
				if line.split('=')[0].strip() in list(FlowInputSettingsDict.keys()):
					exec('self.'+line)
			file.close()

	def complete(self):
		''' Checking if all necessary case attributes are present '''
		from os.path import isfile
		from functions import printOK
		#Only the mesh is absolutely necessary...'
		if isfile(self.MeanFlowFilePath):
			EverythingPresent=True
		else:
			printOK('Set input flow file!')
			EverythingPresent=False
		return EverythingPresent
