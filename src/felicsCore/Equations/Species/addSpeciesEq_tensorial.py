from tensorUtils import *
from ufl import dx
from functions import printWarning, printError, printDebug

def addSpeciesEq(self,fluc,X,mean,species,param):
    '''
    This function builds the weak form of the linearized
    species transport equation in convective form, in the
    tensorial framework.
    '''
    
    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')
        
    if not param.Case.m == 0:
            printWarning('--> Species eq: m > 0 for tensor not validated yet. Treat results with care.')

    # ------------------------ Define the tensorial operators
    J_hat = self._coordinateSystem.J_hat
    
    # ----------------------------------------- Time derivative term
    self.B_vf.add(( fluc.Y(species)*iConj(X)*mean.rho ).ufl_tens*J_hat*dx)
    
    # ----------------------------------------- Advection term
    # This term is integrated by parts
    # Volume term from IbP
    # The volume term seems to introduce a small error (~1e-12) in cartesian nates wrt. previous implementation
    ibp = True
    if ibp:
        self.A_vf.add(( 1j*fluc.Y(species)*iDiv(mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( 1j*mean.Y(species)*iDiv(mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( 1j*mean.Y(species)*iDiv(fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
        # Boundary term from IbP
        self.A_vf.add(( -1j*iDot(self.n,fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(self.n,mean.Y(species)*mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(self.n,mean.Y(species)*fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.all_ds)
    else:
        self.A_vf.add(( -1j*iDot(iGrad(fluc.Y(species)),mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( -1j*iDot(iGrad(mean.Y(species)),mean.rho*fluc.u*iConj(X)) ).ufl_tens*J_hat*dx)
        self.A_vf.add(( -1j*iDot(iGrad(mean.Y(species)),fluc.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*dx)
        


    # ----------------------------------------- Diffusion term
    # NOTE: similarly to what is done in the mom. eq., the diffusion term is integrated by parts but only the 
    # volume part is added to the eqs. --> neglecting the boundary part allows to set a Neumann condition  
    self.A_vf.add(( -1j*mean.D(species)*(iDot(iGrad(fluc.Y(species)),iGrad(iConj(X)))) ).ufl_tens*J_hat*dx)
    self.A_vf.add(( -1j*fluc.D(species)*(iDot(iGrad(mean.Y(species)),iGrad(iConj(X)))) ).ufl_tens*J_hat*dx)
    
    #--------------------------------Reaction
    #reaction = mean.RR_prefactor * mean.rho * (fluc.Y(species) - 2 * fluc.Y(species) * mean.Y(species))\
    #                 + mean.RR_prefactor * fluc.rho * (mean.Y(species) - mean.Y(species) * mean.Y(species))
    ##self.A_vf.add((1j * fluc.omega(species)*iConj(X)).ufl_tens*J_hat*dx)
    #self.A_vf.add((1j * reaction*iConj(X)).ufl_tens*J_hat*dx)
    

    # ----------------------------------------- input output analysis and body forcing
    if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
        self.A_vf.add(( mean.forcing(species)*iConj(X) ).ufl_tens*J_hat*dx)
        

    # ----------------------------------------- BC terms
    for Boundary in param.BCs.getBCsDict()[species]:
        # If forcing is applied at the boundary, and Input-Output mode is on...
        if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
            # First subtract the part added in a few lines above...
            self.A_vf.add(( 1j*iDot(self.n,fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
            # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
            self.A_vf.add(( -1j*iDot(self.n,mean.u*mean.forcing(species)*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
            
        else:
            # Check if boundary condition is Dirichlet or Neumann
            # If it is Dirichlet, the BC is applied in the weak formulation
            if Boundary['type'] in ['Dirichlet']:
                # First subtract the part added in a few lines above...
                self.A_vf.add(( 1j*iDot(self.n, fluc.Y(species)*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                # ... if the value is zero, a treatment is not necessary anyway
                if not Boundary['value'] in [0.0]:
                    printWarning('--> Species eq: Dirichlet BC with non-zero value not validated in tensor framework! Treat results with care.')
                    self.A_vf.add(( -1*iDot(self.n,Boundary['value']*mean.rho*mean.u*iConj(X)) ).ufl_tens*J_hat*self.ds(Boundary['ID']))
                    
            if Boundary['type'] in ['Neumann']:
                if not Boundary['value'] in [0.0]:
                    printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')
