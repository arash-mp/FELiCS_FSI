from ufl import (
				i,
				j,
				k,
				dx,
				Dx,
				inner,
				grad,
				Identity,
				)


def MomentumEqDiffusion(self,fluc,X,mean,param):
	
#------------------------------------------------------------------------------
# VISCOUS TERMS FOR INCOMPRESSIBLE FLOW:
#------------------------------------------------------------------------------
	# I,III: div(mean.mu*grad(fluc.u)) in 2D after Ibp (boundary integral neglected): 
	self.A_imag_vf.add( -self.R*mean.nuTot*Dx(fluc.u[i],j)*Dx(X[i],j)*dx )
	
	#div(fluc.mu*grad(mean.u)) in 2D after Ibp (integral of boundary terms neglected):
	#CAUTION: NOT IMPLEMENTED IN CYLINDRICAL COORDINATES!
	self.A_imag_vf.add( -fluc.nulam*mean.u[i].dx(j)*X[i].dx(j)*dx )
	
#------------------------------------------------------------------------------
# ADDITIONAL TERMS FOR COMPRESSIBLE FLOW:
#------------------------------------------------------------------------------
	if param.Case.HeatTransfer or param.Case.Compressible:
		Iden=Identity( fluc.rho.geometric_dimension() )
		
		# II, III: div(mu*grad^T(fluc.u)) in 2D after Ibp (boundary integral neglected):
		self.A_imag_vf.add( -self.R*mean.nuTot*Dx(fluc.u[j],i)*Dx(X[i],j)*dx )
		
		# IV: div(mean.mu*(div(u)*I)
		self.A_imag_vf.add( mean.nuTot*2.0/3.0*Dx(fluc.u[k],k)*Dx(self.R*X[i],j)*Iden[i,j]*dx ) # Change CA:
		#div(fluc.mu*grad^T(mean.u)) in 2D after Ibp (integral of boundary terms neglected):
		#CAUTION: NOT IMPLEMENTED IN CYLINDRICAL COORDINATES!
		self.A_imag_vf.add( -fluc.nulam*mean.u[j].dx(i)*X[i].dx(j)*dx )
		self.A_imag_vf.add( 2.0/3.0 * fluc.nulam*mean.u[k].dx(k)*Iden[i,j]*X[i].dx(j)*dx )


#------------------------------------------------------------------------------
# ADDITIONAL TERMS FOR CYLINDRICAL COORDINATES, INCOMPRESSIBLE FLOW: 
#------------------------------------------------------------------------------
	if param.Case.CoordinateSystem =='Cylindrical':									#is this necessary? (CA)
		self.A_imag_vf.add(-mean.nuTot/self.R*X[2]*fluc.u[2]*dx)#
		self.A_imag_vf.add(-mean.nuTot.dx(1)*X[2]*fluc.u[2]*dx)
		self.A_imag_vf.add(-mean.nuTot*X[1]/self.R*fluc.u[1]*dx)#
	if not param.Case.m == 0:
		self.A_real_vf.add(-2*mean.nuTot/self.R*self.m*X[2]*fluc.u[1]*dx)#
		self.A_real_vf.add(-X[2]*mean.nuTot.dx(0)*self.m*fluc.u[0]*dx)
		#self.A_imag_vf.add(-X*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)							#III  -> A_real ??? (CA)
		self.A_real_vf.add(-X[2]*mean.nuTot.dx(1)*self.m*fluc.u[1]*dx)
		self.A_imag_vf.add(-mean.nuTot*self.m**2/self.R*inner( X,fluc.u )*dx)#

		self.A_real_vf.add(2*mean.nuTot*self.m/self.R*X[1]*fluc.u[2]*dx)#


			
#------------------------------------------------------------------------------
# ADDITIONAL TERMS FOR CYLINDRICAL COORDINATES, COMPRESSIBLE FLOW (CA): 
#------------------------------------------------------------------------------
		#if True:
		if param.Case.HeatTransfer or param.Case.Compressible:
			
			self.A_imag_vf.add(-mean.nuTot/self.R*fluc.u[1]*X[1]*dx)									#II(1)
			
			
			self.A_imag_vf.add(-2.0/3.0*fluc.u[1]*Dx(mean.nuTot,i)*X[i]*dx)							#IV(1)
			self.A_imag_vf.add(-2.0/3.0*mean.nuTot*Dx(fluc.u[1],i)*X[i]*dx)							#IV(2)
			self.A_imag_vf.add(2.0/3.0*mean.nuTot/self.R*fluc.u[1]*X[1]*dx)							#IV(3)
			
			if not param.Case.m == 0:
				
				self.A_real_vf.add(-self.m*mean.nuTot*Dx(fluc.u[2],i)*X[i]*dx)							#II(2)
				self.A_real_vf.add(self.m/self.R*mean.nuTot*fluc.u[2]*X[1]*dx)							#II(3)

				self.A_real_vf.add(2.0/3.0*self.m*fluc.u[2]*Dx(mean.nuTot,i)*X[i]*dx) 					#IV(4)
				self.A_real_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[2],i)*X[i]*dx)					#IV(5)
				self.A_real_vf.add(-2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[2]*X[1]*dx)				#IV(6)
			
			
			
			
			
			
			
			
#------------------------------------------------------------------------------
#  End of changes in comp. with master branch



	if param.Case.AnalysisMode in ['Input-Output']:
		for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

			# Add physical boundary terms in imaginary part, where the forcing is applied...
			self.A_imag_vf.add(-self.R*mean.nuTot*(-X[0] * (self.n_BC[0] * (fluc.u[0].dx(0)) +
	                                                       		self.n_BC[1] *   (fluc.u[0].dx(1))) -
	                                       		X[1] * (self.n_BC[0] *   (fluc.u[1].dx(0)) +
	                                                       		self.n_BC[1] *   (fluc.u[1].dx(1))) ) *self.ds(boundary_index))
			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_imag_vf.add(-self.R*mean.nuTot*((self.n_BC[0] * (X[0].dx(0)) +
	                   				self.n_BC[1] *   (X[0].dx(1))) * (fluc.u[0]-mean.u_forcing_i[0]) +#jvs
	                  			       (self.n_BC[0] *  (X[1].dx(0)) +
	                   				self.n_BC[1] *   (X[1].dx(1))) * (fluc.u[1]-mean.u_forcing_i[1])) *self.ds(boundary_index)) #jvs

			# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
			self.A_real_vf.add(-self.R*mean.nuTot*((self.n_BC[0] * (X[0].dx(0)) +
	                   				self.n_BC[1] * (X[0].dx(1) )) * (-mean.u_forcing_r[0]) +
	                  			       (self.n_BC[0] * (X[1].dx(0) ) +
	                   				self.n_BC[1] * (X[1].dx(1))) * (-mean.u_forcing_r[1])) *self.ds(boundary_index))


	if param.NumericalScheme in ['Discontinuous Galerkin']:
		                #F_vis_Int = ((self.n_BC[0]('+') * avg(4 / 3 * X[0].dx(0) - 2 / 3 * X[1].dx(1)) +
				#	self.n_BC[1]('+') * avg(X[0].dx(1) + X[1].dx(0))) * jump(fluc['u'][0]) +
				#	(self.n_BC[0]('+') * avg(X[0].dx(1) + X[1].dx(0)) +
				#	 self.n_BC[1]('+') * avg(4 / 3 * X[1].dx(1) - 2 / 3 * X[0].dx(0))) * jump(fluc['u'][1]) -
				#	jump(X[0]) * (n[0]('+') * avg(4 / 3 * fluc['u'][0].dx(0) - 2 / 3 * fluc['u'][1].dx(1)) +
				#	                n[1]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0))) -
				#	jump(Xiu[1]) * (n[0]('+') * avg(fluc['u'][0].dx(1) + fluc['u'][1].dx(0)) +
				#	                n[1]('+') * avg(4 / 3 * fluc['u'][1].dx(1) - 2 / 3 * fluc['u'][0].dx(0))))
				#F_vis_Int *= avg(mean.nuTot)
				F_vis_Int_own = ((   self.n_BC[0]('+') * avg(X[0].dx(0)) +
						     self.n_BC[1]('+') * avg(X[0].dx(1) )) * jump(fluc.u[0]) +
						    (self.n_BC[0]('+') * avg(X[1].dx(0) ) +
						     self.n_BC[1]('+') * avg(X[1].dx(1))) * jump(fluc.u[1]) -
						jump(X[0]) * (self.n_BC[0]('+') * avg(fluc['u'][0].dx(0)) +
						                      self.n_BC[1]('+') * avg(fluc.u[0].dx(1))) -
						jump(X[1]) * (self.n_BC[0]('+') * avg(fluc.u[1].dx(0)) +
						                      self.n_BC[1]('+') * avg(fluc.u[1].dx(1))))
				F_vis_Int_own *= avg(mean.nuTot)
				self.A_imag_vf.add(-F_vis_Int_own *dS)











#####Other stuff
	if  param.Case.TransVelFluc:
		if param.Case.HeatTransfer or param.Case.Compressible:
			self.A_imag_vf.add(self.R*Dx(mean.nuTot,i)*Dx(fluc.u[2],i)*X[2]*dx)						#III(1) (not 100% sure why only in case of dilatation -> needs checking)

			if not param.Case.m == 0:
				self.A_imag_vf.add(-self.m**2/self.R*mean.nuTot*fluc.u[2]*X[2]*dx)					#II(4)
				self.A_imag_vf.add(2/3*self.m**2/self.R*mean.nuTot*fluc.u[2]*X[2]*dx)					#IV(9)

				self.A_real_vf.add(-self.m*mean.nuTot*Dx(fluc.u[i],i)*X[2]*dx)						#II(5)
				self.A_real_vf.add(-self.m/self.R*mean.nuTot*fluc.u[1]*X[2]*dx)						#II(6)
				self.A_real_vf.add(2.0/3.0*self.m*mean.nuTot*Dx(fluc.u[i],i)*X[2]*dx)				#IV(7)
				self.A_real_vf.add(2.0/3.0*self.m/self.R*mean.nuTot*fluc.u[1]*X[2]*dx)				#IV(8)




		if param.Case.AnalysisMode in ['Input-Output']:
			for boundary_index in param.IOResolvent.ForcingBoundaryIndices:

				# Add physical boundary terms in imaginary part, where the forcing is applied...
				self.A_imag_vf.add(-self.R*mean.nuTot*(-X * (self.n_BC[0] *   (fluc.u[2].dx(0)) + self.n_BC[self.ThirdVelCompIndex] *   (fluc.u[2].dx(1)))) *self.ds(boundary_index))

				# Add stabilization terms on imaginary part according to Baumann and Oden JFM 2016 vol 798
				self.A_imag_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (X.dx(0)) +
												self.n_BC[self.ThirdVelCompIndex] * (X.dx(1))) * (fluc.u[2] - mean.ut_forcing_i)) * self.ds(boundary_index))

			# Add stabilization terms on real part according to Baumann and Oden JFM 2016 vol 798
				self.A_real_vf.add(-self.R * mean.nuTot * ((	self.n_BC[0] * (X.dx(0)) +
												self.n_BC[self.ThirdVelCompIndex] * (X.dx(1))) * (-mean.ut_forcing_r)) * self.ds(boundary_index))
