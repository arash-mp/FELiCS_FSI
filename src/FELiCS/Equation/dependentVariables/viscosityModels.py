import numpy as np

def SutherlandFluctuationMean(mean, rho, Ts):
    """
    Calculate the fluctuation viscosity and fluctuation factor using the Sutherland model.

    Parameters:
    mean (object): The mean temperature, density, and laminar viscosity.
    rho (float): The density.
    Ts (float): The fluctuation temperature.

    Returns:
    tuple: A tuple containing the fluctuation viscosity and fluctuation factor.
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
