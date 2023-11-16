from tensor_utils import *
from ufl import (
                i,
                j,
                div,
                inner,
                dx,
                Dx,
                grad,
                Identity,
                conj,
                )

from functions import printWarning, printError

def addSpeciesEq(self,fluc,X,mean,species,param):
    '''
    This function builds the weak form of the linearized
    species transport equation in convective form, in the
    tensorial framework.
    
    NOTE:   The "FLAG_TENS" is used to activate the tensor
            framework in specific terms of the eq. It is possible
            to combine tensorial terms with index-notation ones
            
    NOTE:   How to set "m" in the tensor iGrad and iDiv: the value
            should be the sum of the m for each term in the argument,
            by convention:
                .m = 0 for mean flow quantities
                .m = -m for conj(x)
    '''
    
    # Disclaimers
    if param.NumericalScheme in ['Discontinuous Galerkin']:
        printError('Discontinuous Galerkin not implemented in tensorial framework.')
        
    if not param.Case.m == 0:
            printWarning('--> Species eq: m > 0 for tensor not validated yet. Treat results with care.')


    # ------------------------ Define the tensorial operators
    # Coordinate system
    coord = self.coord_sys
    
    # Vector quantities
    if fluc.u.ufl_shape == 2:
        x_tens = Tensor(as_vector((X[0], X[1], 0.0)), coord)
        u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord)
        u_m = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord)
    
    elif fluc.u.ufl_shape == 3:
        x_tens = Tensor(as_vector((X[0], X[1], X[2])), coord)
        u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord)
        u_m = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord)
        
    else:
        printError('u has neither 2 or 3 dimensions: not implemented.') 
        
    # Scalar quantities
    Y_m = Tensor(mean.Y(species), coord)  
    Y_f = Tensor(fluc.Y(species), coord)
    D_m = Tensor(mean.D(species), coord)
    D_f = Tensor(fluc.D(species), coord)
    rho_m = Tensor(mean.rho, coord)
    rho_f = Tensor(fluc.rho, coord)
    x_tens = Tensor(X, coord)
    if param.Case.AnalysisMode in ['Input-Output']:
        forcing_tens = Tensor(mean.forcing_r(species)+1j*mean.forcing_i(species), coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord) 
    
    
    # ----------------------------------------- Advection term
    # This term is integrated by parts
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # Volume term from IbP
        # The volume term seems to introduce a small error (~1e-12) in cartesian coordinates wrt. previous implementation
        self.A_vf.add(( 1j*Y_f*iDiv(rho_m*u_m*iConj(x_tens),-self.m) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*Y_m*iDiv(rho_m*u_f*iConj(x_tens),0) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( 1j*Y_m*iDiv(rho_f*u_m*iConj(x_tens),0) ).ufl_tens*coord.J_hat*dx)
        # Boundary term from IbP
        self.A_vf.add(( -1j*iDot(nbc_tens,Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(nbc_tens,Y_m*rho_m*u_f*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
        self.A_vf.add(( -1j*iDot(nbc_tens,Y_m*rho_f*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
        
    else:
        # -- > Previous implementation
        I = Identity( fluc.u.geometric_dimension() )
        # Volume term from IbP
        self.A_vf.add(1j * fluc.Y(species)*Dx(mean.u[i]*conj(X)*self.R*mean.rho,j)*I[i,j]*dx)
        self.A_vf.add(1j * mean.Y(species)*Dx(fluc.u[i]*conj(X)*self.R*mean.rho,j)*I[i,j]*dx)
        self.A_vf.add(1j * mean.Y(species)*Dx(mean.u[i]*conj(X)*self.R*fluc.rho,j)*I[i,j]*dx)
        # Boundary term from IbP
        self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*conj(X)*mean.rho*I[i,j]*self.all_ds)
        self.A_vf.add(1j * -self.n_BC[j]*fluc.u[i]*self.R*mean.Y(species)*conj(X)*mean.rho*I[i,j]*self.all_ds)
        self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*mean.Y(species)*conj(X)*fluc.rho*I[i,j]*self.all_ds)
        # Extra term for m > 0 
        if (param.Case.CoordinateSystem in ['Cylindrical']) and (not param.Case.m == 0):
            printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
            self.A_vf += X*param.Case.m*fluc.Y(i_eqn)*mean.ut*dx #possibly conj(X) instead of X

    
    # ----------------------------------------- Time derivative term
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        self.B_vf.add(( Y_f*iConj(x_tens)*rho_m ).ufl_tens*coord.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.B_vf.add( self.R*fluc.Y(species)*conj(X)*mean.rho*dx )    
    

    # ----------------------------------------- Diffusion term
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # NOTE: similarly to what is done in the mom. eq., the diffusion term is integrated by parts but only the 
        # volume part is added to the eqs. --> neglecting the boundary part allows to set a Neumann condition  
        self.A_vf.add(( -1j*D_m*(iDot(iGrad(Y_f,self.m),iGrad(iConj(x_tens),-self.m))) ).ufl_tens*coord.J_hat*dx)
        self.A_vf.add(( -1j*D_f*(iDot(iGrad(Y_m),iGrad(iConj(x_tens),-self.m))) ).ufl_tens*coord.J_hat*dx)
        
    else:
        # -- > Previous implementation
        self.A_vf.add(1j * -mean.D(species)*(self.R*inner(grad(fluc.Y(species)),grad(X)))*dx) # Shouldn't the r be in the derivative???
        self.A_vf.add(1j * -fluc.D(species)*(self.R*inner(grad(mean.Y(species)),grad(X)))*dx)
        # Extra term for m > 0 
        if (param.Case.CoordinateSystem in ['Cylindrical']) and (not param.Case.m == 0):
            printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
            self.A_vf.add(1j * mean.D*X*param.Case.m**2/self.R*fluc.Y(species)*dx)


    # ----------------------------------------- Source terms
    FLAG_TENS = True
    if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
        if FLAG_TENS:
            # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
            self.A_vf.add(( forcing_tens*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
            
        else:
            # -- > Previous implementation
            self.A_vf.add(self.R*mean.forcing_r(species)*conj(X)*dx)
            self.A_vf.add(1j * self.R*mean.forcing_i(species)*conj(X)*dx)
        

    # ----------------------------------------- BC terms
    FLAG_TENS = True
    for Boundary in param.BCs.getBCsDict()[species]:
        # If forcing is applied at the boundary, and Input-Output mode is on...
        if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
            if FLAG_TENS:
                # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
                # First subtract the part added in a few lines above...
                self.A_vf.add(( 1j*iDot(nbc_tens,Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                self.A_vf.add(( -1j*iDot(nbc_tens,u_m*forcing_tens*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                
            else:
                # -- > Previous implementation
                I = Identity( fluc.u.geometric_dimension() )
                # First subtract the part added in a few lines above...
                self.A_vf.add(1j * self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*conj(X)*I[i,j] * self.ds(Boundary['ID']))  ### tlk: Why is therer no density in the equation???
                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_i(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))
                self.A_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_r(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))            
            
        else:
            # Check if boundary condition is Dirichlet or Neumann
            # If it is Dirichlet, the BC is applied in the weak formulation
            if Boundary['type'] in ['Dirichlet']:
                if FLAG_TENS:
                    # First subtract the part added in a few lines above...
                    self.A_vf.add(( 1j*iDot(nbc_tens, Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                    # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                    # ... if the value is zero, a treatment is not necessary anyway
                    if not Boundary['value'] in [0.0]:
                        printWarning('--> Species eq: Dirichlet BC with non-zero value not validated in tensor framework! Treat results with care.')
                        self.A_vf.add(( -1*iDot(nbc_tens,Boundary['value']*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                    
                else:
                    # First subtract the part added in a few lines above...
                    self.A_vf.add(1j * +inner(self.n_BC,mean.u)*self.R*fluc.Y(species)*conj(X)*self.ds(Boundary['ID']))
                    # Then add the respective Boundary condition
                    # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                    # ... if the value is zero, a treatment is not necessary anyway
                    if not Boundary['value'] in [0.0]:
                        self.A_vf.add(-inner(self.n_BC,mean.u)*self.R*Boundary['value']*conj(X)*self.ds(Boundary['ID']))
                    
            if Boundary['type'] in ['Neumann']:
                if not Boundary['value'] in [0.0]:
                    printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')
