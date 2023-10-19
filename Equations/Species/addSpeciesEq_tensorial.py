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

def addSpeciesEq(self,fluc,X,mean,species,param):
    # Define the tensorial operators
    # This should be moved somewhere else
    CoordSys = CoordinateSystem(self.x, param.Case.CoordinateSystem.lower(), mesh_dims = (1, 1, 0))

    # Tensor objects for mean quantities (For vectors we need to pass them as 3D)
    # should tensor objects could be assembled in the meanflow and fluc classes?
    Ymean_T = Tensor(mean.Y(species), CoordSys)
    Dmean_T = Tensor(mean.D(species), CoordSys)
    forcing_T = Tensor(mean.forcing_r(species)+1j*mean.forcing_i(species), CoordSys)
    umean_T = Tensor(as_vector((mean.u[0], mean.u[1], 0.0)), CoordSys)
    rhomean_T = Tensor(mean.rho, CoordSys)

    # Tensor objects for fluctuations
    Yfluc_T = Tensor(fluc.Y(species), CoordSys)
    Dfluc_T = Tensor(fluc.D(species), CoordSys)
    ufluc_T = Tensor(as_vector((fluc.u[0], fluc.u[1], 0.0)), CoordSys)
    rhofluc_T = Tensor(fluc.rho, CoordSys)

    # Tensor object for test function
    X_T = Tensor(X, CoordSys)

    # Surface vector
    nT = Tensor(as_vector((self.n_BC[0], self.n_BC[1], 0.0)), CoordSys)

    # ----------------------------------------- Advection terms
     # -- > Previous implementation
    #I=Identity( fluc.u.geometric_dimension() )
    #self.A_vf.add(1j * fluc.Y(species)*Dx(mean.u[i]*conj(X)*self.R*mean.rho,j)*I[i,j]*dx)
    #self.A_vf.add(1j * mean.Y(species)*Dx(fluc.u[i]*conj(X)*self.R*mean.rho,j)*I[i,j]*dx)
    #self.A_vf.add(1j * mean.Y(species)*Dx(mean.u[i]*conj(X)*self.R*fluc.rho,j)*I[i,j]*dx)
    
    #self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*conj(X)*mean.rho*I[i,j]*self.all_ds)
    #self.A_vf.add(1j * -self.n_BC[j]*fluc.u[i]*self.R*mean.Y(species)*conj(X)*mean.rho*I[i,j]*self.all_ds)
    #self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*mean.Y(species)*conj(X)*fluc.rho*I[i,j]*self.all_ds)
 
    # -- > Tensor implementation (not working for cyl coord)
    self.A_vf.add( ( 1j * Yfluc_T * iDiv(rhomean_T*umean_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * dx)
    self.A_vf.add( (-1j * iDot(nT, Yfluc_T*rhomean_T*umean_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * self.all_ds)

    self.A_vf.add( ( 1j * Ymean_T * iDiv(rhomean_T*ufluc_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * dx)
    self.A_vf.add( (-1j * iDot(nT, Ymean_T*rhomean_T*ufluc_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * self.all_ds)

    self.A_vf.add( ( 1j * Ymean_T * iDiv(rhofluc_T*umean_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * dx)
    self.A_vf.add( (-1j * iDot(nT, Ymean_T*rhofluc_T*umean_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * self.all_ds)
    
    # This should be implemented directly into the tensorial framework
    if param.Case.CoordinateSystem in ['Cylindrical']:
        if not param.Case.m == 0:
            printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
            self.A_vf += X*param.Case.m*fluc.Y(i_eqn)*mean.ut*dx #possibly conj(X) instead of X

    
    # ----------------------------------------- Time derivative terms
    # -- > Previous implementation
    #self.B_vf.add( self.R*fluc.Y(species)*conj(X)*mean.rho*dx )
    # -- > Tensor implementation (not working for cyl coord)
    self.B_vf.add( (Yfluc_T*iConj(X_T)*rhomean_T).ufl_tens * CoordSys.J_hat *dx )
    

    # ----------------------------------------- Diffusion terms
    # Not sure what the term before the weak form is supposed to be
    # -- > Previous implementation
    #self.A_vf.add(1j * -mean.D(species)*(self.R*inner(grad(fluc.Y(species)),grad(X)))*dx) # Shouldn't the r be in the derivative???
    #self.A_vf.add(1j * -fluc.D(species)*(self.R*inner(grad(mean.Y(species)),grad(X)))*dx)
    # -- > Tensor implementation (not working for cyl coord)
    # iInner(a, b) not defined for vector, so replaced by equivalent iDot(a, iConj(b))
    self.A_vf.add( ( -1j* Dmean_T * (iDot(iGrad(Yfluc_T), iGrad(iConj(X_T))))).ufl_tens * CoordSys.J_hat * dx)
    self.A_vf.add( ( -1j* Dfluc_T * (iDot(iGrad(Ymean_T), iGrad(iConj(X_T))))).ufl_tens * CoordSys.J_hat * dx)
    

    # This should be implemented directly into the tensorial framework
    if param.Case.CoordinateSystem in ['Cylindrical']:
        if not param.Case.m == 0:
            printWarning("Cylindrical coordinates for m <> 0 are not validated. Treat results with care")
            self.A_vf.add(1j * mean.D*X*param.Case.m**2/self.R*fluc.Y(species)*dx)
    

    # ----------------------------------------- Source terms
    if param.Case.AnalysisMode == 'Input-Output' and param.IOResolvent.ForcingMode == 'Body':
        # -- > Previous implementation
        #self.A_vf.add(self.R*mean.forcing_r(species)*conj(X)*dx)
        #self.A_vf.add(1j * self.R*mean.forcing_i(species)*conj(X)*dx)
        # -- > Tensor implementation (not working for cyl coord)
        self.A_vf.add( (forcing_T*iConj(X)).ufl_tens * CoordSys.J_hat *dx)
        

    # ----------------------------------------- BC terms
    #I=Identity( fluc.u.geometric_dimension() )
    for Boundary in param.BCs.getBCsDict()[species]:
        # If forcing is applied at the boundary, and Input-Output mode is on...
        if (param.Case.AnalysisMode in ['Input-Output']) and (Boundary['ID'] in param.IOResolvent.ForcingBoundaryIndices) :
            # First subtract the part added in a few lines above...
            # -- > Previous implementation
            #self.A_vf.add(1j * self.n_BC[j]*mean.u[i]*self.R*fluc.Y(species)*conj(X)*I[i,j] * self.ds(Boundary['ID']))  ### tlk: Why is therer no density in the equation???
            # -- > Tensor implementation (not working for cyl coord)
            self.A_vf.add( (1j * iDot(nT, Yfluc_T*rhomean_T*umean_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * self.ds(Boundary['ID']))
            
            # Then add the forcing of the respective species given in the mean flow dict at the respective bounary
            # -- > Previous implementation
            #self.A_vf.add(1j * -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_i(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))
            #self.A_vf.add( -self.n_BC[j]*mean.u[i]*self.R*mean.forcing_r(species)*conj(X)*I[i,j]*self.ds(Boundary['ID']))
            # -- > Tensor implementation (not working for cyl coord)
            self.A_vf.add( (-1 * iDot(nT, umean_T*forcing_T*iConj(X_T)) ).ufl_tens * CoordSys.J_hat * self.ds(Boundary['ID']))
            
        else:
            ## Check if Boundary condition is Dirichlet or Neumann
            if Boundary['type'] in ['Dirichlet']:
                # If it is Dirichlet, the BC is applied in the weak formulation
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
