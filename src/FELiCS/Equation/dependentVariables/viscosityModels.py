import numpy as np

def SutherlandFluctuationMean(mean, rho, Ts):
    """
    Calculate fluctuation viscosity and fluctuation factor using the Sutherland model.

    Parameters
    ----------
    mean : object
        Object containing mean temperature, density, and laminar viscosity (attributes: T, rho, nulam or fieldDict).
    rho : float
        Density value for the fluctuation calculation.
    Ts : float
        Fluctuation temperature.

    Returns
    -------
    nulam : float
        Fluctuation viscosity computed from the Sutherland model.
    fluct : float
        Fluctuation factor for the viscosity.
    """

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
