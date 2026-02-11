from abc import ABC, abstractmethod

class EquationTemplate(ABC):
    """
    A base class for equation templates.

    This abstract base class provides a structure for defining equations 
    within a given coordinate system and with specified parameters. 
    It interfaces with a collection of equations and fluctuation data.

    **Initialize the EquationTemplate object**

    Parameters
    ----------
    index : int
            The index of the equation.
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


    def __init__(self,
    index,
    eqColl,
    fluc,
    X,
    param,
    ):
        """
        Initializes the EquationTemplate instance.

        Parameters
        ----------
        index : int
            The index of the equation.
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuations object.
        X : CoordinateSystem
            The coordinate system object.
        param : dict
            Dictionary of parameters for the equation.
        """

        self.index             = index
        self.fluc              = fluc
        self.X                 = X
        self.param             = param
        self.J_hat             = eqColl._coordinateSystem.J_hat
        self.n                 = eqColl.n
        self.ds                = eqColl.ds
        self.all_ds            = eqColl.all_ds
        self.coordinate_system  = eqColl._coordinateSystem
        self.boundaryHandler   = eqColl.boundaryHandler


    @abstractmethod
    def add_weight_matrix_expression(self):
        """
        Define the weight matrix expression.

        This abstract method must be implemented by subclasses to specify
        how the weight matrix is constructed.
        """
        pass

    @abstractmethod
    def add_linear_expression(self):
        """
        Define the linear expression.

        This abstract method must be implemented by subclasses to specify
        the linear part of the equation.
        """
        pass

    @abstractmethod
    def add_nonlinear_expression(self):
        """
        Define the nonlinear expression.

        This abstract method must be implemented by subclasses to specify
        the nonlinear part of the equation.
        """
        pass

    def add_bilinear_expression(self,
    *args,
    ):
        """
        Optionally define a bilinear expression.

        Parameters
        ----------
        *args : tuple
            Optional arguments required for the bilinear expression.
        """
        pass
