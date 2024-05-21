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

def Epsilon_BilinearTerms(funcSpace, nu, C_mu, Sigma_Epsilon, C1, C2, q, u, k, epsilon, dq, du, dk, depsilon, Xepsilon):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Define nut & Pk & differential form of nut
    '''
    nut = C_mu * k * k / epsilon 
    Pk = 0.5 * nut * Dx(u[i],j) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[i],j) * Dx(u[j], i) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[j], i)    
    nu_e = C_mu / Sigma_Epsilon * k * k / epsilon
    dmudk = C_mu / Sigma_Epsilon * 2 * k * dk / epsilon
    dmude = C_mu / Sigma_Epsilon * k * k / epsilon / epsilon * depsilon * (-1)
    '''
    Term for debugging
    Set depsilon = 0
    '''
    # LNSE_variational =\
    #     1j * depsilon * conj(Xepsilon) * dx
    '''
    Convecting terms and their boundary terms
    '''
    LNSE_variational =\
        - 1j * depsilon * Dx((conj(Xepsilon) * u[i]), i) * dx 
    LNSE_variational +=\
        - 1j * Dx(du[i] * conj(Xepsilon), i) * epsilon * dx 
    LNSE_variational +=\
        + 1j * u[i] * depsilon * n_BC[i] * conj(Xepsilon) * ds 
    LNSE_variational +=\
        + 1j * du[i] * n_BC[i] * epsilon * conj(Xepsilon) * ds 

    '''
    Difussion terms due to viscosity and their boundary terms (commend for zero gradient bc)
    '''
    LNSE_variational +=\
        + 1j * nu * Dx(depsilon, i) * Dx(conj(Xepsilon), i) * dx 
    # LNSE_variational +=\
    #     - 1j * nu * Dx(depsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    '''
    Difussion terms due to turbulent viscosity and their boundary terms (commend for zero gradient bc)
    '''
    LNSE_variational +=\
        + 1j * nu_e * Dx(depsilon, i) * Dx(conj(Xepsilon), i) * dx 
    LNSE_variational +=\
        + 1j * dmudk * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx 
    LNSE_variational +=\
        + 1j * dmude * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx 
    # LNSE_variational +=\
    #     - 1j * nu_e * Dx(depsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    # LNSE_variational +=\
    #     - 1j * dmudk * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    # LNSE_variational +=\
    #     - 1j * dmude * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    '''
    C1 terms(contains Pk, not in weak formulation) 
    '''
    LNSE_variational +=\
        - derivative(1j * C1 * epsilon / k * Pk * conj(Xepsilon) * dx, q, dq) * dq 
    '''
    C2 terms
    '''
    LNSE_variational +=\
        + 1j * C2 * 2 * epsilon / k * depsilon * conj(Xepsilon) * dx \
        - 1j * C2 * epsilon * epsilon / k / k * dk * conj(Xepsilon) * dx 
      
    return LNSE_variational

def Epsilon_LinearTerms(funcSpace,nu, C_mu, Sigma_Epsilon, C1, C2, u, k, epsilon, Xepsilon):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Define nut & Pk & nu_e
    '''
    nut = C_mu * k * k / epsilon 
    Pk = 0.5 * nut * Dx(u[i],j) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[i],j) * Dx(u[j], i) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[i], j) \
        + 0.5 * nut * Dx(u[j], i) * Dx(u[j], i)
    nu_e = C_mu / Sigma_Epsilon * k * k / epsilon
    '''
    Term for debugging
    Set depsilon = 0
    '''
    # NSE_variational =\
    #     0
    '''
    Convecting terms and their boundary terms
    '''
    NSE_variational =\
        - 1j * epsilon * Dx((conj(Xepsilon)* u[i]), i) * dx 
    NSE_variational +=\
        + 1j * u[i] * epsilon * n_BC[i] * conj(Xepsilon) * ds 
    '''
    Difussion terms due to viscosity and their boundary terms (commend for zero gradient bc)
    '''
    NSE_variational +=\
        + 1j * nu * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx 
    # NSE_variational +=\
    #     - 1j * nu * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    '''
    Difussion terms due to turbulent viscosity and their boundary terms (commend for zero gradient bc)
    '''
    NSE_variational +=\
        + 1j * nu_e * Dx(epsilon, i) * Dx(conj(Xepsilon), i) * dx 
    # NSE_variational +=\
    #     - 1j * nu_e * Dx(epsilon, i) * n_BC[i] * conj(Xepsilon) * ds 
    '''
    C1 terms(contains Pk, not in weak formulation) 
    '''
    NSE_variational +=\
        - 1j * C1 *  epsilon / k * Pk * conj(Xepsilon) * dx 
    '''
    C2 terms
    '''
    NSE_variational +=\
        + 1j * C2 * epsilon * epsilon / k * conj(Xepsilon) * dx 
          
    return - NSE_variational
