from ufl import (
    dx
)
from FELiCS.Misc.tensorUtils import (
    iDot,
    iDiv,
    iGrad,
    iConj,
)
from .EquationTemplate      import EquationTemplate
from FELiCS.Misc.logging    import Logger

# Get the logger
logger = Logger.get_logger("felics")

class EnergyPressureEquation(EquationTemplate):
    """
    Class representing the energy-pressure equation in FELiCS.

    This class constructs and manages the weak formulation of the energy 
    conservation equation in a compressible fluid system. The formulation 
    is based on total energy conservation and is expressed in terms of 
    pressure, using the Perfect Gas Law and continuity equations.

    Parameters
    ----------
    eqColl : EquationCollection
        The equation collection object that stores various equations.
    fluc : Fluctuations
        The fluctuations object representing perturbations in the system.
    X : Function
        The solution function representing the unknowns of the equation.
    param : Parameters
        The parameters object containing simulation and physical parameters.

    Raises
    ------
    RuntimeError
        If an unsupported numerical scheme like Discontinuous Galerkin is used.
    """

    def __init__(self, index, eqColl, fluc, X, param):
        """
        Initialize the EnergyPressureEquation object.

        Parameters
        ----------
        eqColl : EquationCollection
            The equation collection object that stores various equations.
        fluc : Fluctuations
            The fluctuations object representing perturbations in the system.
        X : Function
            The solution function representing the unknowns of the equation.
        param : Parameters
            The parameters object containing simulation and physical parameters.

        Raises
        ------
        RuntimeError
            If an unsupported numerical scheme like Discontinuous Galerkin is used.
        """
        # Disclaimer
        if param.Numerics.NumericalScheme in ['Discontinuous Galerkin']:
            logger.error('Discontinuous Galerkin not implemented in tensorial framework.')
            raise Exception('Discontinuous Galerkin not implemented in tensorial framework.')

        # initialize variables in template class
        super().__init__(index, eqColl, fluc, X, param)


    def addWeightMatrixExpression(self, weakForm, mean):
        """
        Add the weight matrix expression to the weak form.

        This function contributes the time derivative terms 
        to the weak formulation.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the expression is added.
        mean : Function
            The mean function representing the average state.
        """
        # Time derivative terms
        # Volume term: -omega*p_f*conj(X)
        weakForm.add((self.fluc.p * iConj(self.X)).ufl_tens * self.J_hat * dx)

    def addNonlinearExpression(self):
        """
        Add the nonlinear expression to the weak form.

        This method is currently not implemented.
        """
        pass

    def addLinearExpression(self, weakForm, mean):
        """
        Add the linear expression to the weak form.

        This function builds the weak form of the linearized energy 
        conservation equation in the tensorial framework. The formulation 
        is based on total energy conservation for a compressible fluid, 
        expressed in terms of pressure by substituting the Perfect Gas Law 
        and the continuity equation. The equation is implemented in 
        primitive variables.

        Parameters
        ----------
        weakForm : Form
            The weak form object to which the expression is added.
        mean : Function
            The mean function representing the average state.
        """
        
        J_hat   = self.J_hat
        X       = self.X
        fluc    = self.fluc

        # ------------------------  Advection terms
        # NOTE: "." denotes the dot product below
        # This includes two groups of terms:
        #   1)  -j*u.grad(p) = -j*(u_f.grad(p_m) + u_m.grad(p_f))
        #   2)  -j*gamma*p*div(u) =  -j*gamma*(p_f*div(u_m) + p_m*div(u_f))
        # both are integrated by parts in the weak form:
        # For 1)
        #   1.vol)  j*[div(conj(X)*u_f)*p_m + j*div(conj(X)*u_m)*p_f]*dx
        #   1.bc)   -j*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
        # Similarly for 2)
        #   2.vol)  j*gamma*[(grad(conj(X)*p_m)).u_f + (grad(conj(X)*p_f)).u_m]*dx
        #   2.bc)  -j*gamma*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
        # NOTE: the 1.bc) and 2.bc) terms are grouped together:
        #   bc)     -j*(1+gamma)*[(conj(X)*u_f*p_m).n_bc + (conj(X)*u_m*p_f).n_bc]*ds
        
        # Add volume integral of pressure gradient term (1)
        weakForm.add(( 1j*iDiv(iConj(X)*mean.u)*fluc.p ).ufl_tens*J_hat*dx)
        weakForm.add(( 1j*iDiv(iConj(X)*fluc.u)*mean.p ).ufl_tens*J_hat*dx)
        # Add volume integral of velocity divergence term (2)
        weakForm.add(( 1j*mean.gamma*iDot(iGrad(iConj(X)*mean.p),fluc.u) ).ufl_tens*J_hat*dx)
        weakForm.add(( 1j*mean.gamma*iDot(iGrad(iConj(X)*fluc.p),mean.u) ).ufl_tens*J_hat*dx)
        # Add boundary integral of (1) and (2)
        weakForm.add(( -1j*(mean.gamma+1)*iDot(mean.u*fluc.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
        weakForm.add(( -1j*(mean.gamma+1)*iDot(fluc.u*mean.p*iConj(X),self.n) ).ufl_tens*J_hat*self.all_ds)
    
        # ------------------------  Thermal diffusion term (Fourier law)
        # NOTE: kappa is the thermal conductivity
        # The term is: -j*(gamma-1)*div(kappa_m*grad(T_f) + kappa_f*grad(T_m))
        # which is integrated by parts,
        # TEMPORARY: use constant Pr definition (This should move to fieldProperty or a handler)
        kappa_m = mean.nuTot*mean.cp/mean.Pr            
        kappa_f = fluc.nulam*mean.cp/mean.Pr
        # Volume term:  
        #   j*(gamma-1)[(grad(conj(X)).(kappa_m*grad(T_f)) + (grad(conj(X)).(kappa_f*grad(T_m))]*dx
        weakForm.add(( -1j*(mean.gamma-1)*iDot(iGrad(iConj(X)),kappa_m*iGrad(fluc.T)) ).ufl_tens*J_hat*dx)
        weakForm.add(( -1j*(mean.gamma-1)*iDot(iGrad(iConj(X)),kappa_f*iGrad(mean.T)) ).ufl_tens*J_hat*dx)
        # Boundary term:
        #   -j*(gamma-1)[(conj(X)*kappa_m*grad(T_f)).n_bc + (conj(X)*kappa_f*grad(T_m)).n_bc]*ds
        #   --> Neglected to impose the proper BC
        # self.A_vf.add(( 1j*(mean.gamma-1)*iDot(iConj(X)*kappa_m*iGrad(fluc.T),self.n) ).ufl_tens*J_hat*self.all_ds)
        # self.A_vf.add(( 1j*(mean.gamma-1)*iDot(iConj(X)*kappa_f*iGrad(mean.T),self.n) ).ufl_tens*J_hat*self.all_ds)
        
        
        # ------------------------  Viscous diffusion term
        # NOTE: kappa is the thermal conductivity
        # The term is: -j*(gamma-1)*(Tau_m*grad(u_f) + Tau_f*grad(u_m))
        # which is integrated by parts,
        # Volume term:
        #   -j*(gamma-1)*[div(Tau) . u*conj(X)]*dx IS it really a minus in front?!?!?!
        weakForm.add(( -1j*(mean.gamma-1)*iDot(iDiv(mean.tau),fluc.u*iConj(X)) ).ufl_tens*J_hat*dx)
        weakForm.add(( -1j*(mean.gamma-1)*iDot(iDiv(fluc.tau),mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
        
        # Boundary term:
        #   +j*(gamma-1)*[u . T*conj(X)]*ds 
        #   --> Neglected to impose the proper BC
        
    

