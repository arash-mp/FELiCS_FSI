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

def Momentum_LinearTerms(funcSpace, nu, C_mu, u, p, k, epsilon, Xu):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Term for debugging
    Set du = 0
    '''
    # NSE_variational =\
    #     0
    '''
    Convecting terms and their boundary terms
    '''
    NSE_variational =\
        - 1j * Dx(conj(Xu[i]) * u[j], j) * u[i] * dx 
    NSE_variational +=\
        + 1j * conj(Xu[i]) * u[j] * u[i] * n_BC[j] * ds 
    '''
    Difussion terms due to viscosity and their boundary terms
    '''
    NSE_variational +=\
        + 1j * nu * Dx(u[i], j) * Dx(conj(Xu[i]), j) * dx \
        + 1j * nu * Dx(u[j], i) * Dx(conj(Xu[i]), j) * dx 
    # NSE_variational +=\
    #     - 1j * nu * Dx(u[i],j) * conj(Xu[j]) * n_BC[i] * ds \
    #     - 1j * nu * Dx(u[j],i) * conj(Xu[i]) * n_BC[j] *ds 
    '''
    Difussion terms due to turbulent viscosity and their boundary terms
    '''
    NSE_variational +=\
        + 1j * C_mu * k * k / epsilon * Dx(u[i], j) * Dx(conj(Xu[i]), j) * dx \
        + 1j * C_mu * k * k / epsilon * Dx(u[j], i) * Dx(conj(Xu[i]), j) * dx 
    # NSE_variational +=\
    #     - 1j * C_mu * k * k / epsilon * Dx(u[i],j) * conj(Xu[j]) * n_BC[i] * ds \
    #     - 1j * C_mu * k * k / epsilon * Dx(u[j],i) * conj(Xu[i]) * n_BC[j] * ds 
    '''
    Pressure term and its boundary term
    '''
    NSE_variational +=\
        - 1j * p * Dx(conj(Xu[i]),i) * dx 
    NSE_variational +=\
        + 1j * p * conj(Xu[i]) * n_BC[i] * ds
    '''
    k term and its boundary term
    '''
    NSE_variational +=\
        - 1j * 2/3 * k * Dx(conj(Xu[i]), i) * dx 
    NSE_variational +=\
        + 1j * 2/3 * k * conj(Xu[i]) * n_BC[i] * ds 
    
     
    

    
                
    return - NSE_variational
