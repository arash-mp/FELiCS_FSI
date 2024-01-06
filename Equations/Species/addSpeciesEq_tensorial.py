from tensorUtils import *
from ufl import dx
from functions import printWarning, printError, printDebug

def addSpeciesEq(self,fluc,X,mean,species,param):
    '''
    This function builds the weak form of the linearized
    species transport equation in convective form, in the
    tensorial framework.
            
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
    #if fluc.u.ufl_shape[0] == 2:
    #    printDebug(param.debug, '--> Species cons. eq: fluctuations and mean flow 2D')
    #    u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), coord)
    #    u_m = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), coord)
    #
    #elif fluc.u.ufl_shape[0] == 3:
    #    printDebug(param.debug, '--> Species cons. eq: fluctuations and mean flow 3D')
    #    u_f = Tensor(as_vector((fluc.u[0], fluc.u[1], fluc.u[2])), coord)
    #    u_m = Tensor(as_vector((mean.u[0], mean.u[1], mean.u[2])), coord)
    #    
    #else:
    #    printError('--> Species-nc eq: u has neither 2 or 3 dimensions: not implemented.') 
    u_f = fluc.u
    u_m = mean.u
        
    # Scalar quantities
    Y_m = mean.Y(species)
    Y_f = fluc.Y(species)
    D_m = mean.D(species)
    D_f = fluc.D(species)
    rho_m = mean.rho
    rho_f = fluc.rho
    #x_tens = Tensor(X, coord, containsTestFunction = True)
    x_tens = X
    if param.Case.AnalysisMode in ['Input-Output']:
        forcing_tens = Tensor(mean.forcing_r(species)+1j*mean.forcing_i(species), coord)
    
    # Boundary normal vector (always 2D)
    nbc_tens = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)),coord)
    
    
    # ----------------------------------------- Advection term
    # This term is integrated by parts
    # Volume term from IbP
    # The volume term seems to introduce a small error (~1e-12) in cartesian coordinates wrt. previous implementation
    self.A_vf.add(( 1j*Y_f*iDiv(rho_m*u_m*iConj(x_tens),self.m) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*Y_m*iDiv(rho_m*u_f*iConj(x_tens),0) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( 1j*Y_m*iDiv(rho_f*u_m*iConj(x_tens),0) ).ufl_tens*coord.J_hat*dx)
    # Boundary term from IbP
    self.A_vf.add(( -1j*iDot(nbc_tens,Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(nbc_tens,Y_m*rho_m*u_f*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)
    self.A_vf.add(( -1j*iDot(nbc_tens,Y_m*rho_f*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.all_ds)

    
    # ----------------------------------------- Time derivative term
    self.B_vf.add(( Y_f*iConj(x_tens)*rho_m ).ufl_tens*coord.J_hat*dx)


    # ----------------------------------------- Diffusion term
    # NOTE: similarly to what is done in the mom. eq., the diffusion term is integrated by parts but only the 
    # volume part is added to the eqs. --> neglecting the boundary part allows to set a Neumann condition  
    self.A_vf.add(( -1j*D_m*(iDot(iGrad(Y_f,self.m),iGrad(iConj(x_tens),self.m))) ).ufl_tens*coord.J_hat*dx)
    self.A_vf.add(( -1j*D_f*(iDot(iGrad(Y_m),iGrad(iConj(x_tens),self.m))) ).ufl_tens*coord.J_hat*dx)


    # ----------------------------------------- Source terms
    if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
        self.A_vf.add(( forcing_tens*iConj(x_tens) ).ufl_tens*coord.J_hat*dx)
        

    # ----------------------------------------- BC terms
    for Boundary in param.BCs.getBCsDict()[species]:
        # If forcing is applied at the boundary, and Input-Output mode is on...
        if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
            # First subtract the part added in a few lines above...
            self.A_vf.add(( 1j*iDot(nbc_tens,Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
            # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
            self.A_vf.add(( -1j*iDot(nbc_tens,u_m*forcing_tens*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
            
        else:
            # Check if boundary condition is Dirichlet or Neumann
            # If it is Dirichlet, the BC is applied in the weak formulation
            if Boundary['type'] in ['Dirichlet']:
                # First subtract the part added in a few lines above...
                self.A_vf.add(( 1j*iDot(nbc_tens, Y_f*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                # If BC value is 0, then the weak formulation throws an error, therefore check if it is zero...
                # ... if the value is zero, a treatment is not necessary anyway
                if not Boundary['value'] in [0.0]:
                    printWarning('--> Species eq: Dirichlet BC with non-zero value not validated in tensor framework! Treat results with care.')
                    self.A_vf.add(( -1*iDot(nbc_tens,Boundary['value']*rho_m*u_m*iConj(x_tens)) ).ufl_tens*coord.J_hat*self.ds(Boundary['ID']))
                    
            if Boundary['type'] in ['Neumann']:
                if not Boundary['value'] in [0.0]:
                    printError('So far only homogeneous Neumann conditions are implemented... Please either change to another BC or - even better -  implement it yourself and upload your well documented implementation to gitlab...')
