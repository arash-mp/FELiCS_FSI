from ufl import (
				lhs,
				rhs,
				)

class WeakForm():
#	def __init__(self):
#		self.weakForm=0
	def add(self, input_vf):
		if hasattr(self,'__weakForm__'):
			self.__weakForm__ += input_vf
		else:
			self.__weakForm__ = input_vf
	def subtract(self,input_vf):
		if hasattr(self,'__weakForm__'):
			self.__weakForm__ -= input_vf
		else:
			self.__weakForm__ = -input_vf

	def lhsIsZero(self):
		if self.isZero():
			return True
		else:
			temp=lhs(self.__weakForm__)
			if len(temp.arguments())<2:
				return True
			else:
				return False

	def rhsIsZero(self):
		if self.isZero():
			return True
		else:
			temp=rhs(self.__weakForm__)
			#print(temp)
			#print(temp.arguments())
			#input(len(temp.arguments()))
			if len(temp.arguments())<1:
				return True
			else:
				return False

	def isZero(self):
		if hasattr(self,'__weakForm__'):
			return False
		else:
			return True
	@property
	def lhs(self):
		return lhs(self.__weakForm__)
	@property
	def rhs(self):
		return rhs(self.__weakForm__)

	@property
	def weakForm(self):
		return self.__weakForm__


	def analyseWeakForm(self):

		if hasattr(self,'__weakForm__'):
			print('###### Weak form:')
			print(self.weakForm)
			print('###### Weak form arguments:')
			print(self.weakForm.arguments())
			print('###### Number of weak form arguments:')
			print(len(self.weakForm.arguments()))
		else:
			print('Weak form is zero')
