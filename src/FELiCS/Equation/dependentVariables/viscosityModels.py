#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |  
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |  
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
import numpy as np

def sutherland_fluctuation_mean(mean,
                                rho,
                                Ts,
                                ):
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
        mean_T = mean.field_dict['T']
        mean_rho = mean.field_dict['rho']
        mean_nulam = mean.field_dict['nulam']
    else:
        mean_T = mean.T
        mean_rho = mean.rho
        mean_nulam = mean.nulam

    fluct = (mean_T + 3 * Ts) / (2 * (mean_T + Ts)) * (-1 * rho / mean_rho)
    nulam = mean_nulam * fluct
    return nulam, fluct
