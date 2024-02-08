
import numpy as np

def SutherlandFluctuationMean(mean, rho,Ts):
    #if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != rho.shape[0]:
    if False:
        local_T = mean.fieldDict['T']
        local_rho = mean.fieldDict['rho']
        local_nulam = mean.fieldDict['nulam']
    else:
        local_T = mean.T
        local_rho = mean.rho
        local_nulam = mean.nulam
    fluct = (local_T + 3 * Ts) / (2 * (local_T + Ts)) * (-1 * rho / local_rho)
    nulam = local_nulam * fluct
    return nulam, fluct
