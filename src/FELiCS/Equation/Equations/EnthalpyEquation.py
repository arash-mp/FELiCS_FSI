from ufl import (
    dx,
    conj,
    Identity,
    i,
    j,
    k,
    Dx,
    as_tensor,
    inner,
    grad,
    dot,
    outer,
    transpose,
    Constant,
)
from FELiCS.Misc.tensorUtils import (
    Tensor,
    as_vector,
    iInner,
    iDot,
    iDiv,
    iGrad,
    iConj,
    iOuter,
    iT,
    iIdentity,
)




class EnthalpyEquation(EquationTemplate):
    """Class representing the enthalpy conservation equation."""

    def __init__(self, eqColl, fluc, X, param):
        """
        Initialize the EnthalpyEquation object.

        Parameters:
        -----------
        eqColl : EquationCollection
            The equation collection object.
        fluc : Fluctuations
            The fluctuations object.
        X : Function
            The function representing the mesh coordinates.
        param : Parameters
            The parameters object.
        """
        # Disclaimer
        if param.NumericalScheme in ['Discontinuous Galerkin']:
            printError('Discontinuous Galerkin not implemented in tensorial framework.')
    
        # initialize variables in template class
        super().__init__(eqColl, fluc, X, param)


    def addWeightMatrixExpression(self, weakForm, mean):
        """
        Add the weight matrix expression to the weak form.

        Parameters:
        -----------
        weakForm : Form
            The weak form object.
        mean : MeanFields
            The mean fields object.
        """
        #  Time derivative terms
        weakForm.add((mean.rho * self.fluc.h * iConj(self.X)).ufl_tens * self.J_hat * dx)
        weakForm.add((-1 * self.fluc.p * iConj(self.X)).ufl_tens * self.J_hat * dx)
    
     
    def addNonlinearExpression(self):
        """
        Add the nonlinear expression to the weak form.
        """
        pass

    def addLinearExpression(self, weakForm, mean):
        """
        Add the linear expression to the weak form.

        Parameters:
        -----------
        weakForm : Form
            The weak form object.
        mean : MeanFields
            The mean fields object.
        """
        '''
        This function builds the weak form of the linearized
        enthalpy conservation equation, in tensorial framework.
        '''
        printDebug(self.param.debug, "Adding transport equation for enthalpy in all mesh internal directions")

        J_hat = self.J_hat
        X = self.X
        fluc = self.fluc
        
            
        # ------------------------  Advection terms
        # Add volume integral of advection terms that remain after partial integration:
        weakForm.add((1j * iDiv(iConj(X) * mean.rho * mean.u) * fluc.h).ufl_tens * J_hat * dx)
        weakForm.add((1j * iDiv(iConj(X) * fluc.rho * mean.u) * mean.he).ufl_tens * J_hat * dx)
        weakForm.add((1j * iDiv(iConj(X) * mean.rho * fluc.u) * mean.he).ufl_tens * J_hat * dx)
        # Add boundary integrals resulting from said partial integration:
        weakForm.add((-1j * iDot(mean.rho * mean.u * fluc.h * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
        weakForm.add((-1j * iDot(fluc.rho * mean.u * mean.he * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
        weakForm.add((-1j * iDot(mean.rho * fluc.u * mean.he * iConj(X), self.n)).ufl_tens * J_hat * self.all_ds)
            
        # ------------------------  Diffusion terms
        weakForm.add((-1j * iDot(mean.alpha * iGrad(fluc.h), iGrad(iConj(X)))).ufl_tens * J_hat * dx)
        weakForm.add((-1j * iDot(fluc.alpha * iGrad(mean.he), iGrad(iConj(X)))).ufl_tens * J_hat * dx)
