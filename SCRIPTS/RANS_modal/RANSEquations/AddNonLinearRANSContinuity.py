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
                conj
                )

def Continuity_BilinearTerms(funcSpace, du, dp, Xp):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Term for debugging
    Set dp = 0
    '''
    # LNSE_variational =\
    #     1j * dp * conj(Xp) * dx
    LNSE_variational =\
        - 1j * du[i] * Dx(conj(Xp),i) * dx\
        + 1j * du[i] * n_BC[i] * conj(Xp) * ds
        
    return LNSE_variational

def Continuity_LinearTerms(funcSpace, u, Xp):
    n_BC=FacetNormal(funcSpace.mesh)
    '''
    Term for debugging
    Set dp = 0
    '''
    # NSE_variational =\
    #     0
    NSE_variational =\
        - 1j * u[i] * Dx(conj(Xp),i) * dx\
        + 1j * u[i] * n_BC[i] * conj(Xp) *ds
            
    return - NSE_variational
