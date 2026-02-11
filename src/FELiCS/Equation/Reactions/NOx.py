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

import pathlib
from FELiCS.Misc.functions import loadCSV

from dolfinx.fem import (
                        Function
                        )
from ufl import (
                dx,
)
from scipy import interpolate
import numpy as np
class NOx:
    """
    Class for computing NO and NO2 source terms based on pre-tabulated data.

    This class uses a CSV table to interpolate reaction source terms for NO and NO2
    as a function of the mixture fraction or equivalence ratio (`phi`). These interpolated
    values are used to add source terms to the weak form of a finite element formulation.

    **Initialize the NOx object**

    Parameters
    ----------
    P : any
        Placeholder for configuration or context parameter (currently unused).
    """


    def __init__(self,
    P,
    ):
        """
        Initializes the NOx instance by loading tabulated NOx source term data.

        Parameters
        ----------
        P : any
            Placeholder parameter (not currently used).
        """
        tablePath=str(pathlib.Path(__file__).parent.absolute())+'/NOxTable.csv'
        self.__Table=loadCSV(tablePath)



    def add_source_to_weak_form(self,
    weakform,
    dQ_threshold=None,
    ):
        """
        Add NO and NO2 source terms to the weak form based on interpolated data.

        Interpolates pre-tabulated NO and NO2 production rates as a function of
        the local `phi` field and assigns them as source terms in the weak form.
        Optionally filters source terms based on a heat release threshold.

        Parameters
        ----------
        weakform : object
            The weak form structure, expected to include test functions `v_NO`, `v_NO2`,
            and solution fields `phi` and `T`.
        dQ_threshold : float, optional
            A threshold for the mean heat release rate; source terms are suppressed
            where the temperature `T` is below 1000 K.

        Returns
        -------
        eq : ufl.Form
            The assembled weak form with NO and NO2 source terms included.
        """
        interpolationNO = interpolate.interp1d(
        self.__Table['phi'],
        self.__Table['omega_NO_pf'],
        )
        interpolationNO2 = interpolate.interp1d(
        self.__Table['phi'],
        self.__Table['omega_NO2_pf'],
        )
        omegaNO=Function(weakform.phi.function_space)
        print(np.max(weakform.phi.vector[:]))
        omegaNO.vector()[:]=interpolationNO(weakform.phi.vector[:])
        omegaNO2=Function(weakform.phi.function_space())
        omegaNO2.vector()[:]=interpolationNO2(weakform.phi.vector[:])
        print(len(omegaNO2.vector[:]))
        print(len(weakform.T.vector[:]))
        for i in range(len(omegaNO2.vector()[:])):
            if weakform.T.vector[i]<1000:
                omegaNO.vector[i]=0
                omegaNO2.vector[i]=0

        eq =  weakform.v_NO * omegaNO * dx
        eq += weakform.v_NO2 * omegaNO2 * dx
        return eq
