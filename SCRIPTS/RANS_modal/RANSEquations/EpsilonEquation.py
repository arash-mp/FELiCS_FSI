from ufl import (
                dx, 
                Dx,
                i,
                j,
                conj,
                derivative,
                )


from FELiCS.Equation.Equations.EquationTemplate import EquationTemplate

class EpsilonEquation(EquationTemplate):

    def __init__(self, equation, X_epsilon, fluc_u, fluc_k, fluc_epsilon):

        # initialize variables in template class
        super().__init__(equation,0,0,0)

        self.n_BC         = equation.n_BC 

        self.X_epsilon    = X_epsilon

        self.fluc_u       = fluc_u
        self.fluc_k       = fluc_k
        self.fluc_epsilon = fluc_epsilon


    def addWeightMatrixExpression(self,weakForm,mean):
        # Time derivative term
        weakForm.add( self.fluc_epsilon * conj(self.X_epsilon) *dx)

    def addNonlinearExpression(self):
        pass

    def addLinearExpression(self,weakForm,mean):
        '''
        This function builds the weak form of the linearized
        momentum conservation equation, in tensorial framework.
        The current version is based on the addMomentumEq.py function
        that was previously defined using the index notation of ufl.
        For now all terms from the momentum equation have been grouped 
        into a single file instead of being in serparate files.
        '''
        n_BC          = self.n_BC 
        ds            = self.all_ds

        Xepsilon      = self.X_epsilon

        u             = mean._fieldDict['u']
        k             = mean._fieldDict['k']
        epsilon       = mean._fieldDict['epsilon']
        nu            = mean._fieldDict['nulam']
        C_mu          = mean.C_mu
        Sigma_Epsilon = mean.sigma_epsilon
        C1            = mean.C_1epsilon
        C2            = mean.C_2epsilon

        du            = self.fluc_u
        dk            = self.fluc_k
        depsilon      = self.fluc_epsilon


        '''
        Define nut & Pk & differential form of nut
        '''
        nut   = C_mu * k * k / epsilon 
        Pk    = 0.5 * nut * Dx(u[i],j) * Dx(u[i],j) \
              + 0.5 * nut * Dx(u[i],j) * Dx(u[j],i) \
              + 0.5 * nut * Dx(u[j],i) * Dx(u[i],j) \
              + 0.5 * nut * Dx(u[j],i) * Dx(u[j],i)    
        nu_e  = C_mu / Sigma_Epsilon * k * k / epsilon
        dmudk = C_mu / Sigma_Epsilon * 2 * k * dk / epsilon
        dmude = C_mu / Sigma_Epsilon * k * k / epsilon / epsilon * depsilon * (-1)

        nut_fluc  = 2.* C_mu * k * dk / epsilon  + C_mu * k * k * depsilon / epsilon /epsilon
        Pk_fluc   = 0.5 * nut_fluc * (Dx(u[i],j) * Dx(u[i],j) + Dx(u[i],j) * Dx(u[j],i) + Dx(u[j],i) * Dx(u[i],j) + Dx(u[j],i) * Dx(u[j],i))
        Pk_fluc  += 0.5 * nut      * (Dx(du[i],j)* Dx(u[i],j) + Dx(du[i],j)* Dx(u[j],i) + Dx(du[j],i)* Dx(u[i],j) + Dx(du[j],i)* Dx(u[j],i))
        Pk_fluc  += 0.5 * nut      * (Dx(u[i],j) * Dx(du[i],j)+ Dx(u[i],j) * Dx(du[j],i)+ Dx(u[j],i) * Dx(du[i],j)+ Dx(u[j],i) * Dx(du[j],i))
        '''
        Convecting terms and their boundary terms
        '''
        weakForm.add(\
            - 1j * depsilon * Dx((conj(Xepsilon) * u[i]), i) * dx )
        weakForm.add(\
            - 1j * Dx(du[i] * conj(Xepsilon), i) * epsilon * dx )
        weakForm.add(\
            + 1j * u[i] * depsilon * n_BC[i] * conj(Xepsilon) * ds )
        weakForm.add(\
            + 1j * du[i] * n_BC[i] * epsilon * conj(Xepsilon) * ds )
    
        '''
        Difussion terms due to viscosity and their boundary terms (commend for zero gradient bc)
        '''
        weakForm.add(\
            + 1j * nu * Dx(depsilon, i) * Dx(conj(Xepsilon), i) * dx )
        #TODO: remove this term if Neumann is set for depsilon 
        weakForm.add(\
    	     - 1j * nu * Dx(depsilon, i) * n_BC[i] * conj(Xepsilon) * ds) 
        '''
        Difussion terms due to turbulent viscosity and their boundary terms (commend for zero gradient bc)
        '''
        weakForm.add(\
            + 1j * nu_e * Dx(depsilon, i) * Dx(conj(Xepsilon), i) * dx )
        weakForm.add(\
            + 1j * dmudk * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx )
        weakForm.add(\
            + 1j * dmude * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx )
        #TODO: remove this term if Neumann is set for depsilon 
        weakForm.add(\
            - 1j * nu_e * Dx(depsilon, i) * n_BC[i] * conj(Xepsilon) * ds )
        weakForm.add(\
            - 1j * dmudk * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds )
        weakForm.add(\
            - 1j * dmude * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds ) 
        '''
        C1 terms(contains Pk, not in weak formulation) 
        '''
        #weakForm.add(\
        #    - derivative(1j * C1 * epsilon / k * Pk * conj(Xepsilon) * dx, q, dq) * dq )
        weakForm.add(\
             (- 1j * C1 * (depsilon / k + epsilon * dk / k / k) * Pk \
              - 1j * C1 * epsilon / k * Pk_fluc) * conj(Xepsilon) * dx)
        '''
        C2 terms
        '''
        weakForm.add(\
            + 1j * C2 * 2 * epsilon / k * depsilon * conj(Xepsilon) * dx \
            - 1j * C2 * epsilon * epsilon / k / k * dk * conj(Xepsilon) * dx )
          


