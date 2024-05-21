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
                conj,
                derivative
                )

def K_BilinearTerms(funcSpace,nu, C_mu, Sigma_k, q, dq, u, k, epsilon, du, dk, depsilon, Xk):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Define nut & Pk 
    '''
    nut = C_mu * k * k / epsilon 
    Pk = 0.5 * nut * Dx(u[i],j) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[i],j) * Dx(u[j], i) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[j], i)
    '''
    Term for debugging
    Set dk = 0
    '''
    # LNSE_variational =\
    #     1j * dk * conj(Xk) * dx
    '''
    Convecting terms and their boundary terms
    '''
    LNSE_variational =\
        - 1j * Dx(du[i] * conj(Xk), i) * k * dx 
    LNSE_variational +=\
        - 1j * Dx(u[i] * conj(Xk), i) * dk * dx
    LNSE_variational +=\
        + 1j * du[i] * n_BC[i] * k * conj(Xk) * ds
    LNSE_variational +=\
        + 1j * u[i] * n_BC[i] * dk * conj(Xk) * ds 
    '''
    Difussion terms due to viscosity and their boundary terms (commend for zero gradient bc)
    '''
    LNSE_variational +=\
        + 1j * nu * Dx(dk, i) * Dx(conj(Xk), i) * dx 
    # LNSE_variational +=\
    #     - 1j * nu * Dx(dk, i) * n_BC[i] * conj(Xk) * ds     
    '''
    Difussion terms due to turbulent viscosity and their boundary terms (commend for zero gradient bc)
    '''
    LNSE_variational +=\
        + 1j * C_mu / Sigma_k * k * k / epsilon * Dx(dk, i) * Dx(conj(Xk), i) * dx \
        + 1j * C_mu / Sigma_k * 2 * k * dk / epsilon * Dx(k, i) * Dx(conj(Xk), i) * dx 
    LNSE_variational +=\
        - 1j * C_mu / Sigma_k * k * k / epsilon / epsilon * depsilon * Dx(k, i) * Dx(conj(Xk), i) * dx 
    # LNSE_variational +=\
    #     - 1j * C_mu / Sigma_k * k * k / epsilon * Dx(dk, i) * n_BC[i] * conj(Xk) * ds  \
    #     - 1j * C_mu / Sigma_k * 2 * k * dk / epsilon * Dx(k, i) * n_BC[i] * conj(Xk) * ds  
    # LNSE_variational +=\
    #     + 1j * C_mu / Sigma_k * k * k / epsilon / epsilon * depsilon * Dx(k, i) * n_BC[i] * conj(Xk) * ds  
    '''
    Pk terms (not in weak form)
    '''
    LNSE_variational +=\
        - derivative(1j * Pk * conj(Xk) * dx, q, dq) * dq
    '''
    epsilon term
    '''
    LNSE_variational +=\
        + 1j * depsilon * conj(Xk) * dx 
            
    return LNSE_variational

