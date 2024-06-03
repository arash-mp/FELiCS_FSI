from ufl import (
                dx, 
                Dx,
                i,
                j,
                conj,
                derivative,
                )


from FELiCS.Equation.Equations.EquationTemplate import EquationTemplate

class KEquation(EquationTemplate):

    def __init__(self, equation, X_k, fluc_u, fluc_k, fluc_epsilon):

        # initialize variables in template class
        super().__init__(equation,0,0,0)

        self.n_BC         = equation.n_BC 

        self.X_k          = X_k

        self.fluc_u       = fluc_u
        self.fluc_k       = fluc_k
        self.fluc_epsilon = fluc_epsilon


    def addWeightMatrixExpression(self,weakForm,mean):
        # Time derivative term
        weakForm.add( self.fluc_k * conj(self.X_k) *dx)

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
        n_BC     = self.n_BC 
        ds       = self.all_ds

        Xk       = self.X_k

        u        = mean._fieldDict['u']
        k        = mean._fieldDict['k']
        epsilon  = mean._fieldDict['epsilon']
        nu       = mean._fieldDict['nulam']
        C_mu     = mean.C_mu
        Sigma_k  = mean.sigma_k

        du       = self.fluc_u
        dk       = self.fluc_k
        depsilon = self.fluc_epsilon

        '''
        Define nut & Pk 
        '''
        nut = C_mu * k * k / epsilon 
        Pk = 0.5 * nut * Dx(u[i],j) * Dx(u[i],j) \
           + 0.5 * nut * Dx(u[i],j) * Dx(u[j],i) \
           + 0.5 * nut * Dx(u[j],i) * Dx(u[i],j) \
           + 0.5 * nut * Dx(u[j],i) * Dx(u[j],i)
        nut_fluc  = 2.* C_mu * k * dk / epsilon  + C_mu * k * k * depsilon / epsilon /epsilon
        Pk_fluc   = 0.5 * nut_fluc * (Dx(u[i],j) * Dx(u[i],j) + Dx(u[i],j) * Dx(u[j],i) + Dx(u[j],i) * Dx(u[i],j) + Dx(u[j],i) * Dx(u[j],i))
        Pk_fluc  += 0.5 * nut      * (Dx(du[i],j)* Dx(u[i],j) + Dx(du[i],j)* Dx(u[j],i) + Dx(du[j],i)* Dx(u[i],j) + Dx(du[j],i)* Dx(u[j],i))
        Pk_fluc  += 0.5 * nut      * (Dx(u[i],j) * Dx(du[i],j)+ Dx(u[i],j) * Dx(du[j],i)+ Dx(u[j],i) * Dx(du[i],j)+ Dx(u[j],i) * Dx(du[j],i))
        '''
        Convecting terms and their boundary terms
        '''
        weakForm.add(\
            - 1j * Dx(du[i] * conj(Xk), i) * k * dx) 
        weakForm.add(\
            - 1j * Dx(u[i] * conj(Xk), i) * dk * dx)
        weakForm.add(\
            + 1j * du[i] * n_BC[i] * k * conj(Xk) * ds)
        weakForm.add(\
            + 1j * u[i] * n_BC[i] * dk * conj(Xk) * ds)
        '''
        Difussion terms due to viscosity and their boundary terms (commend for zero gradient bc)
        '''
        weakForm.add(\
            + 1j * nu * Dx(dk, i) * Dx(conj(Xk), i) * dx )
        weakForm.add(\
             - 1j * nu * Dx(dk, i) * n_BC[i] * conj(Xk) * ds)     
        '''
        Difussion terms due to turbulent viscosity and their boundary terms (commend for zero gradient bc)
        '''
        weakForm.add(\
            + 1j * C_mu / Sigma_k * k * k / epsilon * Dx(dk, i) * Dx(conj(Xk), i) * dx \
            + 1j * C_mu / Sigma_k * 2 * k * dk / epsilon * Dx(k, i) * Dx(conj(Xk), i) * dx )
        weakForm.add(\
            - 1j * C_mu / Sigma_k * k * k / epsilon / epsilon * depsilon * Dx(k, i) * Dx(conj(Xk), i) * dx )
        #TODO: remove the first line of this term if Neumann is set for dk & depsilon
        weakForm.add(\
            - 1j * C_mu / Sigma_k * k * k / epsilon * Dx(dk, i) * n_BC[i] * conj(Xk) * ds  \
            - 1j * C_mu / Sigma_k * 2 * k * dk / epsilon * Dx(k, i) * n_BC[i] * conj(Xk) * ds)
        weakForm.add(\
            + 1j * C_mu / Sigma_k * k * k / epsilon / epsilon * depsilon * Dx(k, i) * n_BC[i] * conj(Xk) * ds)
        '''
        Pk terms (not in weak form)
        '''
        #weakForm.add(\
        #    - derivative(1j * Pk * conj(Xk) * dx, q, dq) * dq)
        weakForm.add(\
             - 1j * Pk_fluc * conj(Xk) * dx)
        '''
        epsilon term
        '''
        weakForm.add(\
            + 1j * depsilon * conj(Xk) * dx )
            

