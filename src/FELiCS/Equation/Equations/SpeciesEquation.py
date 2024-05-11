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


from FELiCS.Misc.functions import printWarning, printError, printDebug

from .EquationTemplate import EquationTemplate

class SpeciesEquation(EquationTemplate):

    def __init__(self,eqColl,fluc,X,species,param):
 
        # Disclaimer
        if param.NumericalScheme in ['Discontinuous Galerkin']:
            printError('Discontinuous Galerkin not implemented in tensorial framework.')
    
        # initialize variables in template class
        super().__init__(eqColl,fluc,X,param)

        self.species = species


    def addWeightMatrixExpression(self,weakForm,mean):
        # Time derivative term
        weakForm.add(( self.fluc.Y(self.species)*iConj(self.X)*mean.rho ).ufl_tens*self.J_hat*dx)

    def addNonlinearExpression(self):
        pass

    def addLinearExpression(self,weakForm,mean):
        '''
        This function builds the weak form of the linearized
        species transport equation in convective form, in the
        tensorial framework.
        '''
       
        param   = self.param
        J_hat   = self.J_hat
        species = self.species
        X       = self.X
        fluc    = self.fluc

            
        if not param.Case.m == 0:
                printWarning('--> Species eq: m > 0 for tensor not validated yet. Treat results with care.')
        
        
        # ----------------------------------------- Advection term
        # This term is integrated by parts
        # Volume term from IbP
        # The volume term seems to introduce a small error (~1e-12) in cartesian nates wrt. previous implementation
        ibp = True
        if ibp:
            weakForm.add(( 1j*fluc.Y(species)*iDiv(mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( 1j*mean.Y(species)*iDiv(mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( 1j*mean.Y(species)*iDiv(fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
            # Boundary term from IbP
            weakForm.add(( -1j*iDot(self.n,fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
            weakForm.add(( -1j*iDot(self.n,mean.Y(species)*mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
            weakForm.add(( -1j*iDot(self.n,mean.Y(species)*fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
        else:
            weakForm.add(( -1j*iDot(iGrad(fluc.Y(species)),mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( -1j*iDot(iGrad(mean.Y(species)),mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*dx)
            weakForm.add(( -1j*iDot(iGrad(mean.Y(species)),fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
            
    
    
        # ----------------------------------------- Diffusion term
        # NOTE: similarly to what is done in the mom. eq., the diffusion term is integrated by parts but only the 
        # volume part is added to the eqs. --> neglecting the boundary part allows to set a Neumann condition  
        weakForm.add(( -1j*mean.D(species)*(iDot(iGrad(fluc.Y(species)),iGrad(iConj(X)))) ).ufl_tens*J_hat*dx)
        weakForm.add(( -1j*fluc.D(species)*(iDot(iGrad(mean.Y(species)),iGrad(iConj(X)))) ).ufl_tens*J_hat*dx)
        
        #--------------------------------Reaction
        #reaction = mean.RR_prefactor * mean.rho * (fluc.Y(species) - 2 * fluc.Y(species) * mean.Y(species))\
        #                 + mean.RR_prefactor * fluc.rho * (mean.Y(species) - mean.Y(species) * mean.Y(species))
        ##weakForm.add((1j * fluc.omega(species)*iConj(X)).ufl_tens*J_hat*dx)
        #weakForm.add((1j * reaction*iConj(X)).ufl_tens*J_hat*dx)


        # ----------------------------------------- input output analysis and body forcing
        if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
            weakForm.add(( mean.forcing(species)*iConj(X) ).ufl_tens*J_hat*dx)
            
    
        # ----------------------------------------- BC terms
        for Boundary in param.BCs.getBCsDict()[species]:
            # If forcing is applied at the boundary, and Input-Output mode is on...
            if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
                # First subtract the part added in a few lines above...
                weakForm.add(( 1j*iDot(self.n,fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                weakForm.add(( -1j*iDot(self.n,mean.u*mean.forcing(species)*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                
            else:
                # Check if boundary condition is Dirichlet or Neumann
                # If it is Dirichlet, the BC is applied in the weak formulation
                if Boundary['type'] in ['Dirichlet']:
                    # First subtract the part added in a few lines above...
                    weakForm.add(( 1j*iDot(self.n, fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                    # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                    # ... if the value is zero, a treatment is not necessary anyway
                    if not Boundary['value'] in [0.0]:
                        printWarning('--> Species eq: Dirichlet BC with non-zero value not validated in tensor framework! Treat results with care.')
                        weakForm.add(( -1*iDot(self.n,Boundary['value']*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                        
                if Boundary['type'] in ['Neumann']:
                    if not Boundary['value'] in [0.0]:
                        printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')
    
