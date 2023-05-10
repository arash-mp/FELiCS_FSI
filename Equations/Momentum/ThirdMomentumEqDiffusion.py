from ufl import (
				i,
				j,
				k,
				dx,
				Dx,
				grad,
				inner,
				)
def ThirdMomentumEqDiffusion(self,fluc,X,mean,param):


	#self.A_imag_vf.add(-self.R*mean.nuTot*inner( grad(X) , grad(fluc.ut) )*dx)				#I
	

	if param.Case.CoordinateSystem =='Cylindrical':									#is this necessary? (CA)
		self.A_imag_vf.add(-mean.nuTot/self.R*X*fluc.ut*dx)									#I
		self.A_imag_vf.add(-mean.nuTot.dx(1)*X*fluc.ut*dx)									#III
	if not param.Case.m == 0:
		self.A_imag_vf.add(-mean.nuTot*self.m**2/self.R*X*fluc.ut*dx)						#I
		self.A_real_vf.add(-2*mean.nuTot/self.R*self.m*X*fluc.u[1]*dx)						#I
		self.A_real_vf.add(-X*mean.nuTot.dx(0)*self.m*fluc.u[0]*dx)							#III
		#self.A_imag_vf.add(-X*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)							#III  -> A_real ??? (CA)
		self.A_real_vf.add(-X*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)
		
	#if True:
	if param.Case.HeatTransfer or param.Case.Compressible:
		self.A_imag_vf.add(self.R*Dx(mean.nuTot,i)*Dx(fluc.ut,i)*X*dx)						#III(1) (not 100% sure why only in case of dilatation -> needs checking)
		
		if not param.Case.m == 0: 
			self.A_imag_vf.add(-self.m**2/self.R*mean.nuTot*fluc.ut*X*dx)					#II(4)
			self.A_imag_vf.add(2/3*self.m**2/self.R*mean.nuTot*fluc.ut*X*dx)					#IV(9)
		
			self.A_real_vf.add(-self.m*mean.nuTot*Dx(fluc.u[i],i)*X*dx)						#II(5)
			self.A_real_vf.add(-self.m/self.R*mean.nuTot*fluc.u[1]*X*dx)						#II(6)
			self.A_real_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[i],i)*X*dx)				#IV(7)
			self.A_real_vf.add(2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[1]*X*dx)				#IV(8)


	

	if param.Case.AnalysisMode in ['Input-Output']:
		for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

			# Add physical boundary terms in imaginary part, where the forcing is applied...
			self.A_imag_vf.add(-self.R*mean.nuTot*(-X * (self.n_BC[0] *   (fluc.ut.dx(0)) + self.n_BC[self.ThirdVelCompIndex] *   (fluc.ut.dx(1)))) *self.ds(boundary_index))

			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_imag_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (X.dx(0)) +
											self.n_BC[self.ThirdVelCompIndex] * (X.dx(1))) * (fluc.ut - mean.ut_forcing_i)) * self.ds(boundary_index))

		# Add stabilization terms on real part according to Baumann and Oden JFM 2016 vol 798
			self.A_real_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (X.dx(0)) +
											self.n_BC[self.ThirdVelCompIndex] * (X.dx(1))) * (-mean.ut_forcing_r)) * self.ds(boundary_index))


