from ufl import (
                i,
                j,
                div,
                inner,
                dx,
                Dx,
                Identity,
                grad,
                )
from tensor_utils import *

def addSpeciesConservativeEq(self,fluc,X,mean,species,param):
    
    '''
    This function builds the weak form of the linearized
    species conservation equation (conservative form), in 
    the tensorial framework.
    
    NOTE:   The "FLAG_TENS" is used to activate the tensor
            framework in specific terms of the eq. It is possible
            to combine tensorial terms with index-notation ones
            
    THIS EQUATION IS NOT YET IMPLEMENTED IN TENSOR FRAMEWORK BECAUSE WE
    NEED A VALIDATION CASE FOR IT!
    '''
    
    from functions import printDebug, printError, printWarning
    
    ## ------------------------ Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED "mesh_dims" and "as_vector(x1, x2, 0.0)" to cartesian coords
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord_sys)  # HARDCODED FOR 2D perturbations so far
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)  # HARDCODED FOR 2D perturbations so far
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # CYL COORD IN FELiCS DEFINED AS (Z, R, PHI)!
        printWarning("Cylindrical coordinates for tensor species-cons not validated yet. Treat results with care.")
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord_sys)
        
    else:
        printError('Coord. syst not yet implemented in tensor framework.')
    
    rhoY_mean_tens = Tensor(mean.rhoY(species), coord_sys)
    D_mean_tens = Tensor(mean.D(species), coord_sys)
    rhoY_fluc_tens = Tensor(fluc.rhoY(species), coord_sys)
    D_fluc_tens = Tensor(fluc.D(species), coord_sys)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys) 
    rho_mean_tens = Tensor(mean.rho, coord_sys)
    rho_fluc_tens = Tensor(fluc.rho, coord_sys)
    x_tens = Tensor(X, coord_sys)
    
    ## ------------------------  Advection terms
    # addSpeciesConservativeEqAdvection(self,fluc,X,mean,species,param)
    
    # This term is integrated by parts
    FLAG_TENS = False
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # Volume term from IbP
        # The volume term seems to introduce a small error (~1e-12) in cartesian coordinates wrt. previous implementation
        printWarning('Tensor form debugging')
        self.A_vf.add( ( 1j*iDot(rhoY_fluc_tens*u_mean_tens, iGrad(iConj(x_tens), -self.m)) ).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add( ( 1j*iDot(rhoY_mean_tens*u_fluc_tens, iGrad(iConj(x_tens), -self.m)) ).ufl_tens*coord_sys.J_hat*dx)
        
        # Boundary term from IbP
        self.A_imag_vf.add(-self.n_BC[j]*fluc.u[i]*self.R*mean.Y(species)*X*mean.rho*I[i,j]*self.all_ds)
        self.A_imag_vf.add(-self.n_BC[j]*mean.u[i]*self.R*fluc.rhoY(species)*X*I[i,j]*self.all_ds)
        
    else:
        # -- > Previous implementation 
        I=Identity( fluc.u.geometric_dimension() )

        self.A_imag_vf.add( mean.u[i]*fluc.rhoY(species)*Dx(X*self.R,j)*I[i,j]*dx)
        self.A_imag_vf.add( fluc.u[i]*mean.Y(species)*mean.rho*Dx(X*self.R,j)*I[i,j]*dx)

        self.A_imag_vf.add(-self.n_BC[j]*fluc.u[i]*self.R*mean.Y(species)*X*mean.rho*I[i,j]*self.all_ds)
        self.A_imag_vf.add(-self.n_BC[j]*mean.u[i]*self.R*fluc.rhoY(species)*X*I[i,j]*self.all_ds)

        if param.Case.CoordinateSystem in ['Cylindrical']:
            if not param.Case.m == 0:
                printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
                self.A_real_vf += X*param.Case.m*fluc.Y(i_eqn)*mean.ut*dx
 
    
    ## ------------------------  Time derivative terms
    # addSpeciesConservativeEqTimeDerivative(self,fluc,X,mean,species)
 
    self.B_real_vf.add( self.R*fluc.rhoY(species)*X*dx )
 
    
    ## ------------------------  Diffusion terms
    # addSpeciesConservativeEqDiffusion(self,fluc,X,mean,species,param)
 
    self.A_imag_vf.add(-mean.D(species)*(self.R*inner(grad(fluc.Y(species)),grad(X)))*dx) # Shouldn't the r be in the derivative???
    self.A_imag_vf.add(-fluc.D(species)*(self.R*inner(grad(mean.Y(species)),grad(X)))*dx)

    if param.Case.CoordinateSystem in ['Cylindrical']:
        if not param.Case.m == 0:
            printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
            self.A_imag_vf.add(mean.D*X*param.Case.m**2/self.R*fluc.Y(species)*dx)
    
    
    ## ------------------------  Source terms
    # addSpeciesConservativeEqSourceTerms(self,X,mean,species,param)
    
    if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
        self.A_real_vf.add(self.R*mean.forcing_r(species)*X*dx)
        self.A_imag_vf.add(self.R*mean.forcing_i(species)*X*dx)
  
  
    ## ------------------------  BC terms
    # addSpeciesConservativeEqBCs(self,fluc,X,mean,species,param)
    
    I=Identity( fluc.u.geometric_dimension() )
    for Boundary in param.BCs.getBCsDict()[species]:
        # If forcing is applied at the boundary, and Input-Output mode is on...
        if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
            # First subtract the part added in a few lines above...
            self.A_imag_vf.add(self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*X*I[i,j]*self.ds(Boundary['ID']))  ### tlk: Why is therer no density in the equation???
            # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
            self.A_imag_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_i(species)*X*I[i,j]*self.ds(Boundary['ID']))
            self.A_real_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_r(species)*X*I[i,j]*self.ds(Boundary['ID']))
        else:
            ## Check if Boundary condition is Dirichlet or Neumann
            if Boundary['type'] in ['Dirichlet']:
                # If it is Dirichlet, the BC is applied in the weak formulation
                # First subtract the part added in a few lines above...
                self.A_imag_vf.add(+inner(self.n_BC,mean.u)*self.R*fluc.Y(species)*X*self.ds(Boundary['ID']))
                # Then add the respective Boundary condition
                # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                # ... if the value is zero, a treatment is not necessary anyway
                if not Boundary['value'] in [0.0]:
                    self.A_real_vf.add(-inner(self.n_BC,mean.u)*self.R*Boundary['value']*X*self.ds(Boundary['ID']))
            if Boundary['type'] in ['Neumann']:
                if not Boundary['value'] in [0.0]:
                    printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')

