from .SettingsClass import Settings
class NumericsSettingsClass(Settings):
	def __init__(self):
		'''Initializing the numerics settings class
		'''
		super().__init__()
		self._settingsKind = 'Numerics'
		SettingsDict=self.getAllSettingsDict()
		for key in list(SettingsDict.keys()):
		
			if type(SettingsDict[key]['default']) in [dict]:
				quotation=''
			else:
				quotation='\"'
			tempStr='self.'+key+'='+quotation+str(SettingsDict[key]['default'])+quotation
			exec(tempStr)

	def getAllSettingsDict(self):
		'''Function returning all numerics settings with default values'''
		SettingsDict={\
			'nSolut':{'datatype':int,'default':3},
			'EigenValueGuess':{'datatype':list,'default':[1.0]},
			'nCPU':{'datatype':int,'default':1},
			'Schemes':{'datatype':dict,'default':{'u':'CG'},'options':['CG']},
			'PolynomialOrder':{'datatype':dict,'default':{'u':'2'},'options':[1,2]},
			'LinearAlgebraSolver':{'datatype':str,'default':'python'},
			'Preconditioner':{'datatype':str,'default':'None'},
		}
		return SettingsDict

	def importSettings(self,settingFilePath):
		''' Loading Mean Flow parameters from file '''
		SettingsDict=self.getAllSettingsDict()
		if not settingFilePath =='':
			file = open(settingFilePath)

			#Read whole file
			lines = file.readlines()
			# Add every line of the file as an attribute to the object
			
			for line in lines:
				beforeEqualSign = line.split('=')[0].strip()
				afterEqualSign = line.split('=')[1].strip()
				if beforeEqualSign in list(SettingsDict.keys()):
					exec('self.tempVariable = ' + afterEqualSign)
					if type(SettingsDict[beforeEqualSign]['default']) == type(self.tempVariable):
						exec('self.'+line)
					elif SettingsDict[beforeEqualSign]['datatype'] == list and type(self.tempVariable) in [float,int,complex]:
						self.tempVariable = list([self.tempVariable])
						exec('self.'+beforeEqualSign + '=' + str(self.tempVariable))
			del self.tempVariable
			file.close()

	def complete(self,variableList):
		''' Checking if all necessary case attributes are present '''
		#Check inputs for completeness and correctness ...'
		from functions import printWarning
		EverythingPresent=True
		CaseSettingsDict=self.getAllSettingsDict()
		# Check if schemes are chosen correctly for every variable
		for variable in variableList:
			if not variable in list(self.Schemes.keys()) or \
			   not self.Schemes[variable] in CaseSettingsDict['Schemes']['options']:
				printWarning("Schemes are not correctly chosen...")
				EverythingPresent=False
				break
		for variable in variableList:
			if not variable in list(self.PolynomialOrder.keys()) or \
			   not self.PolynomialOrder[variable] in CaseSettingsDict['PolynomialOrder']['options']:
				printWarning("Polynomial orders are not correctly chosen...")
				EverythingPresent=False
				break

		return EverythingPresent
