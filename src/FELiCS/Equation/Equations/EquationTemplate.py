from abc import ABC, abstractmethod

class EquationTemplate(ABC):

    def __init__(self,eqColl,fluc,X,param):
 
        # variables    
        self.J_hat  = eqColl._coordinateSystem.J_hat
        self.param  = param
        self.fluc   = fluc
        self.X      = X
        self.n      = eqColl.n
        self.all_ds = eqColl.all_ds
        self.ds     = eqColl.ds


    @abstractmethod
    def addWeightMatrixExpression(self):
        pass


    @abstractmethod
    def addLinearExpression(self):
        pass

 
    @abstractmethod
    def addNonlinearExpression(self):
        pass
