'''
The functions to add nonlinear momentum fuction to WeakForm class.
It is the basic function to calculate mean field.
!!! So far it is not completed yet.
-- Xiuyang take response to develop it.
'''
from ufl import (
                dx, 
                ds,
                FacetNormal,
                Dx,
                i,
                j,
                conj
                )

def Momentum_BilinearTerms(funcSpace,nu, C_mu, u, k, epsilon, du, dp, dk, depsilon, Xu):
    n_BC=FacetNormal(funcSpace.mesh)

    '''
    Difussion terms due to turbulent viscosity and their boundary terms
    '''
    LNSE_variational +=\
        + 1j * C_mu * k * k / epsilon * Dx(du[i],j) * Dx(conj(Xu[i]), j) * dx \
        + 1j * C_mu * k * k / epsilon * Dx(du[j],i) * Dx(conj(Xu[i]), j) * dx \
    # LNSE_variational +=\
    #     - 1j * C_mu * k * k / epsilon * Dx(du[i],j) * conj(Xu[i]) * n_BC[j] * ds \
    #     - 1j * C_mu * k * k / epsilon * Dx(du[j],i) * conj(Xu[i]) * n_BC[j] * ds \
    LNSE_variational +=\
        + 2j * C_mu * k / epsilon * dk * Dx(u[i],j) * Dx(conj(Xu[i]), j) * dx \
        + 2j * C_mu * k / epsilon * dk * Dx(u[j],i) * Dx(conj(Xu[i]), j) * dx \
    # LNSE_variational +=\
    #     - 2j * C_mu * k / epsilon * dk * Dx(u[i],j) * conj(Xu[j]) * n_BC[i] * ds \
    #     - 2j * C_mu * k / epsilon * dk * Dx(u[j],i) * conj(Xu[i]) * n_BC[j] * ds \
    LNSE_variational +=\
        - 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[i],j) * Dx(conj(Xu[i]), j) * dx \
        - 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[j],i) * Dx(conj(Xu[i]), j) * dx \
    # LNSE_variational +=\
    #     + 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[i],j) * conj(Xu[j]) * n_BC[i] * ds \
    #     + 1j * C_mu * k * k / epsilon / epsilon * depsilon * Dx(u[j],i) * conj(Xu[i]) * n_BC[j] * ds \
    '''
    k term and its boundary term
    '''
    LNSE_variational +=\
        - 1j * 2/3 * dk * Dx(conj(Xu[i]), i) *dx 
    LNSE_variational +=\
        + 1j * 2/3 * dk * conj(Xu[i]) * n_BC[i] * ds 
        

    return LNSE_variational
