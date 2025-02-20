from abc import ABC, abstractmethod

class EquationTemplate(ABC):
    """
    A base class for equation templates.

    Parameters
    ----------
    eqColl : EquationCollection
        The equation collection object.
    fluc : Fluctuations
        The fluctuations object.
    X : CoordinateSystem
        The coordinate system object.
    param : dict
        The parameters for the equation.

    Attributes
    ----------
    J_hat : ndarray
        The J_hat attribute of the coordinate system.
    param : dict
        The parameters for the equation.
    fluc : Fluctuations
        The fluctuations object.
    X : CoordinateSystem
        The coordinate system object.
    n : int
        The value of n.
    all_ds : list
        The list of all ds.
    ds : list
        The list of ds.

    """

    def __init__(self, eqColl, fluc, X, param):
        self.J_hat = eqColl._coordinateSystem.J_hat
        self.param = param
        self.fluc = fluc
        self.X = X
        self.n = eqColl.n
        self.all_ds = eqColl.all_ds
        self.ds = eqColl.ds
        self.coordinateSystem  = eqColl._coordinateSystem


    @abstractmethod
    def addWeightMatrixExpression(self):
        """
        Abstract method to add the weight matrix expression.
        """
        pass

    @abstractmethod
    def addLinearExpression(self):
        """
        Abstract method to add the linear expression.
        """
        pass

    @abstractmethod
    def addNonlinearExpression(self):
        """
        Abstract method to add the nonlinear expression.
        """
        pass
