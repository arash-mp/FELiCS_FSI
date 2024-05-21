from ufl import (
                dx, 
                Dx,
                i,
                j,
                conj,
                )


from FELiCS.Equation.Equations.EquationTemplate import EquationTemplate

class MomentumEquation(EquationTemplate):

    def __init__(self, equation, X_u, fluc_u, fluc_k, fluc_epsilon):

        # initialize variables in template class
        super().__init__(equation,0,0,0)

        self.n_BC         = equation.n_BC 

        self.Xu           = X_u

        self.fluc_u       = fluc_u
        self.fluc_k       = fluc_k
        self.fluc_epsilon = fluc_epsilon


    def addWeightMatrixExpression(self,weakForm,mean):
        # Time derivative term
        pass

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

        Xu       = self.Xu

        u        = mean._fieldDict['u']
        k        = mean._fieldDict['k']
        epsilon  = mean._fieldDict['epsilon']
        C_mu     = mean.C_mu

        du       = self.fluc_u
        dk       = self.fluc_k
        depsilon = self.fluc_epsilon


        '''
        Difussion terms due to turbulent viscosity and their boundary terms
        '''
        weakForm.add(\
            + 1j * C_mu * k * k / epsilon * Dx(du[i],j) * Dx(conj(Xu[i]), j) * dx \
            + 1j * C_mu * k * k / epsilon * Dx(du[j],i) * Dx(conj(Xu[i]), j) * dx) 
        weakForm.add(\
            + 2j * C_mu * k / epsilon * dk * Dx(u[i],j) * Dx(conj(Xu[i]), j) * dx \
            + 2j * C_mu * k / epsilon * dk * Dx(u[j],i) * Dx(conj(Xu[i]), j) * dx) 
        weakForm.add(\
            - 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[i],j) * Dx(conj(Xu[i]), j) * dx \
            - 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[j],i) * Dx(conj(Xu[i]), j) * dx) 
        '''
        k term and its boundary term
        '''
        weakForm.add(\
            - 1j * 2/3 * dk * Dx(conj(Xu[i]), i) *dx )
        weakForm.add(\
            + 1j * 2/3 * dk * conj(Xu[i]) * n_BC[i] * ds )
        

