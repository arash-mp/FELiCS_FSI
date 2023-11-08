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
    # Define the tensorial operators
    if param.Case.CoordinateSystem =='Cartesian':
        # HARDCODED "mesh_dims" and "as_vector(x1, x2, 0.0)" to cartesian coords
        coord_sys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord_sys)
    
    elif param.Case.CoordinateSystem =='Cylindrical':
        # CYL COORD IN FELiCS DEFINED AS (Z, R, PHI)!
        printWarning("Cylindrical coordinates for tensor species-cons not validated yet. Treat results with care.")
        coord_sys = CoordinateSystem(self.x, "cylindricalfelics", mesh_dims = (1, 1, 0))
        u_fluc_tens = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord_sys)
        u_mean_tens = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord_sys)
        
    else:
        printError('Coord. syst not yet implemented in tensor framework.')
        
    if not param.Case.m == 0:
            printWarning('m > 0 for tensor species-cons not validated yet. Treat results with care.')
        
    Y_mean_tens = Tensor(mean.Y(species), coord_sys)
    forcing_tens = Tensor(mean.forcing_r(species)+1j*mean.forcing_i(species), coord_sys)
    D_mean_tens = Tensor(mean.D(species), coord_sys)
    Y_fluc_tens = Tensor(fluc.Y(species), coord_sys)
    D_fluc_tens = Tensor(fluc.D(species), coord_sys)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), coord_sys)
    rho_mean_tens = Tensor(mean.rho, coord_sys)
    rho_fluc_tens = Tensor(fluc.rho, coord_sys)
    x_tens = Tensor(X, coord_sys)

    # ----------------------------------------- Advection term
    # This term is integrated by parts
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # Volume term from IbP
        self.A_vf.add( ( 1j * Y_fluc_tens * iDiv(rho_mean_tens*u_mean_tens*iConj(x_tens), -self.m)).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add( ( 1j * Y_mean_tens * iDiv(rho_mean_tens*u_fluc_tens*iConj(x_tens), 0)).ufl_tens*coord_sys.J_hat*dx)
        self.A_vf.add( ( 1j * Y_mean_tens * iDiv(rho_fluc_tens*u_mean_tens*iConj(x_tens), 0)).ufl_tens*coord_sys.J_hat*dx)
        
        # Boundary term from IbP
        self.A_vf.add( (-1j * iDot(nbc_tens, Y_fluc_tens*rho_mean_tens*u_mean_tens*iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.all_ds)
        self.A_vf.add( (-1j * iDot(nbc_tens, Y_mean_tens*rho_mean_tens*u_fluc_tens*iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.all_ds)
        self.A_vf.add( (-1j * iDot(nbc_tens, Y_mean_tens*rho_fluc_tens*u_mean_tens*iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.all_ds)
        
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
        self.B_vf.add( (Y_fluc_tens*iConj(x_tens)*rho_mean_tens).ufl_tens * coord_sys.J_hat *dx )
        
    else:
        # -- > Previous implementation
        self.B_vf.add( self.R*fluc.Y(species)*conj(X)*mean.rho*dx )    
    

    # ----------------------------------------- Diffusion term
    FLAG_TENS = True
    if FLAG_TENS:
        # -- > Tensor implementation (NOT TESTED FOR CYL. COORDS!)
        # NOTE: Not sure what the term before the weak form is supposed to be
        # Shouldn't there be a rho_mean here? Or is it absorbed into the D?
        self.A_vf.add( ( -1j* D_mean_tens * (iDot(iGrad(Y_fluc_tens, self.m), iGrad(iConj(x_tens), -self.m)))).ufl_tens * coord_sys.J_hat * dx)
        self.A_vf.add( ( -1j* D_fluc_tens * (iDot(iGrad(Y_mean_tens), iGrad(iConj(x_tens), -self.m)))).ufl_tens * coord_sys.J_hat * dx)
        
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
            printWarning("Input-Output analysis with Body forcing not validated yet in tensor framework. Treat results with care")
            self.A_vf.add( (forcing_tens*iConj(x_tens)).ufl_tens * coord_sys.J_hat *dx)
            
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
                printWarning("Input-Output analysis with boundary forcing not validated yet in tensor framework. Treat results with care")
                # First subtract the part added in a few lines above...
                self.A_vf.add( (1j * iDot(nbc_tens, Y_fluc_tens*rho_mean_tens*u_mean_tens*iConj(x_tens)) ).ufl_tens * coord_sys.J_hat * self.ds(Boundary['ID']))
                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                self.A_vf.add( (-1 * iDot(nbc_tens, u_mean_tens*forcing_tens*iConj(x_tens)) ).ufl_tens * coord_sys.J_hat * self.ds(Boundary['ID']))
                
            else:
                # -- > Previous implementation
                I = Identity( fluc.u.geometric_dimension() )
                # First subtract the part added in a few lines above...
                self.A_vf.add(1j * self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*conj(X)*I[i,j] * self.ds(Boundary['ID']))  ### tlk: Why is therer no density in the equation???
                # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
                self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_i(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))
                self.A_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_r(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))            
            
        else:
            # Check if Boundary condition is Dirichlet or Neumann
            # If it is Dirichlet, the BC is applied in the weak formulation
            if Boundary['type'] in ['Dirichlet']:
                if FLAG_TENS:
                    printWarning('Dirichlet BC implementation in species conservation eq. not validated in tensor framework. Treat results with care.')
                    
                    # First subtract the part added in a few lines above...
                    self.A_vf.add((1j * iDot(nbc_tens, Y_fluc_tens*rho_mean_tens*u_mean_tens*iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.ds(Boundary['ID']))
                    # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                    # ... if the value is zero, a treatment is not necessary anyway
                    if not Boundary['value'] in [0.0]:
                        self.A_vf.add((1j * iDot(nbc_tens, Boundary['value']*rho_mean_tens*u_mean_tens*iConj(x_tens))).ufl_tens*coord_sys.J_hat*self.ds(Boundary['ID']))
                    
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
