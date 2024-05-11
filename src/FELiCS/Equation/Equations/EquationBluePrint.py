from abc import ABC, abstractmethod

class EquationBluePrint(ABC):

    @abstractmethod
    def addWeightMatrixExpression(self):
        pass


    @abstractmethod
    def addLinearExpression(self):
        pass

 
    @abstractmethod
    def addNonlinearExpression(self):
        pass
