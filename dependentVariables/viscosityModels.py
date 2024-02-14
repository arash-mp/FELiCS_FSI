
import numpy as np

def SutherlandFluctuationMean(mean, rho,Ts):
    #if isinstance(mean.T, np.ndarray) and mean.T.shape[0] != rho.shape[0]:
    if False:
        mean_T = mean.fieldDict['T']
        mean_rho = mean.fieldDict['rho']
        mean_nulam = mean.fieldDict['nulam']
    else:
        mean_T = mean.T
        mean_rho = mean.rho
        mean_nulam = mean.nulam
    fluct = (mean_T + 3 * Ts) / (2 * (mean_T + Ts)) * (-1 * rho / mean_rho)
    nulam = mean_nulam * fluct
    return nulam, fluct
