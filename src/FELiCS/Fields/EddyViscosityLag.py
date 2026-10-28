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
"""
Time-lagged response of the frozen eddy viscosity (Case.EddyViscosityLag).

With TurbulenceModel "File" the coherent Reynolds stress follows the strain of
the perturbation instantly, tau~ = -2 mu_t S(u').  This module lets it respond
with a delay.  Two one-parameter models, both linear and both exact in the
frequency domain for perturbations ~ exp(-i lambda t):

"Relaxation"   T d(tau~)/dt + tau~ = -2 mu_t S(u')      ->  mu_t / (1 - i lambda T)
"Delay"        tau~(t) = -2 mu_t S(u'(t - T))           ->  mu_t exp(i lambda T)

The relaxation lags the strain by arctan(omega T) (at most 90 deg) and reduces
the amplitude by cos of that angle; the delay lags by omega T without
attenuation.  The factor multiplies the eddy viscosity only where it acts on
the perturbation: the perturbation stress of the momentum equation and the
perturbation wall traction of the rigid-body force.  The mean flow, the base
residual and the moving base state are unchanged, and TimeScale 0 recovers
the frozen model exactly.

The factor contains the eigenvalue, so the eigenproblem is nonlinear in
lambda: the user gives Eigenvalue (the value used in the factor) and repeats
the run with the eigenvalue found until it no longer changes (typically 2-4
runs; see documentation/GoverningEquations/eddy_viscosity_lag.md and
tools/eddy_viscosity_lag_iterate.py).

The time scale can be restricted to a region (Region): "All", "Wake" (x
beyond TrailingEdgeX), or "Band" (XRange[0] < x < XRange[1]), with smooth
ramps of half width RampWidth.  Outside the region the eddy viscosity stays
frozen.  2-D Cartesian cases only.
"""

# Standard library
import json

# Third-party
import numpy as np
import ufl
from dolfinx import fem

# Local
from FELiCS.Misc.logging import Logger, log_and_raise

logger = Logger.get_logger("felics")

MODELS = ("None", "Relaxation", "Delay")
REGIONS = ("All", "Wake", "Band")
REQUIRED_KEYS = {
    "All": ("Model", "TimeScale", "Eigenvalue", "Region"),
    "Wake": ("Model", "TimeScale", "Eigenvalue", "Region", "TrailingEdgeX", "RampWidth"),
    "Band": ("Model", "TimeScale", "Eigenvalue", "Region", "XRange", "RampWidth"),
}


def smoothstep(values):
    """Cubic smoothstep of the clipped argument (numpy)."""
    clipped = np.clip(values, 0.0, 1.0)
    return clipped * clipped * (3.0 - 2.0 * clipped)


class EddyViscosityLag:
    """
    Lagged response of the frozen eddy viscosity to the perturbation.

    Parameters
    ----------
    param : object
        Global parameter object; reads ``Case.EddyViscosityLag``.
    mean_flow : MeanFlowClass
        Mean flow after its viscosities are initialised (needs the mesh and
        the eddy-viscosity field).

    Methods
    -------
    factor()
        UFL expression (complex scalar) multiplying the eddy viscosity in the
        perturbation stress and traction; 1 for the frozen model.
    """

    def __init__(
        self,
        param,
        mean_flow,
    ):
        settings = dict(getattr(param.Case, "EddyViscosityLag", None) or {"Model": "None"})
        self._model = settings.get("Model", "None")
        self._mean_flow = mean_flow
        self._time_scale_field = None
        self._time_scale = 0.0
        self._eigenvalue = 0.0
        self._region = "All"
        if self._model not in MODELS:
            log_and_raise(
                logger,
                f"Case.EddyViscosityLag.Model must be one of {list(MODELS)}, got '{self._model}'.",
                ValueError,
            )
        if self._model == "None":
            unknown = set(settings) - {"Model"}
            if unknown:
                log_and_raise(
                    logger,
                    f"Case.EddyViscosityLag: Model 'None' takes no other keys, got {sorted(unknown)}.",
                    ValueError,
                )
            return
        self._check_settings(settings, param)
        self._time_scale = float(settings["TimeScale"])
        real, imag = settings["Eigenvalue"]
        self._eigenvalue = complex(float(real), float(imag))
        self._region = settings["Region"]
        if self._region != "All":
            self._time_scale_field = self._masked_time_scale(settings)
        logger.info(
            f"Eddy-viscosity lag: model {self._model}, time scale {self._time_scale}, "
            f"eigenvalue {self._eigenvalue} in the factor, region {self._region}, "
            f"settings {json.dumps(settings)}."
        )

    # ------------------------------------------------------------------
    # protected
    # ------------------------------------------------------------------
    def _check_settings(
        self,
        settings,
        param,
    ):
        region = settings.get("Region")
        if region not in REGIONS:
            log_and_raise(
                logger,
                f"Case.EddyViscosityLag.Region must be one of {list(REGIONS)}, got '{region}'.",
                ValueError,
            )
        required = REQUIRED_KEYS[region]
        missing = [key for key in required if key not in settings]
        unknown = [key for key in settings if key not in required]
        if missing or unknown:
            log_and_raise(
                logger,
                f"Case.EddyViscosityLag (Region {region}) needs exactly the keys "
                f"{list(required)}; missing {missing}, unknown {unknown}.",
                ValueError,
            )
        if float(settings["TimeScale"]) < 0.0:
            log_and_raise(logger, "Case.EddyViscosityLag.TimeScale must be >= 0.", ValueError)
        eigenvalue = settings["Eigenvalue"]
        if not (isinstance(eigenvalue, (list, tuple)) and len(eigenvalue) == 2):
            log_and_raise(
                logger,
                "Case.EddyViscosityLag.Eigenvalue must be [real, imaginary] (FELiCS "
                "convention, perturbations ~ exp(-i lambda t)).",
                ValueError,
            )
        if region != "All" and float(settings["RampWidth"]) <= 0.0:
            log_and_raise(logger, "Case.EddyViscosityLag.RampWidth must be positive.", ValueError)
        if region == "Band":
            x_range = settings["XRange"]
            if not (isinstance(x_range, (list, tuple)) and len(x_range) == 2 and x_range[0] < x_range[1]):
                log_and_raise(
                    logger,
                    "Case.EddyViscosityLag.XRange must be [x_start, x_end] with x_start < x_end.",
                    ValueError,
                )
        if param.Case.CoordinateSystem != "Cartesian" or int(param.Case.nDim) != 2:
            log_and_raise(
                logger,
                "Case.EddyViscosityLag is implemented for 2-D Cartesian cases only.",
                NotImplementedError,
            )

    def _masked_time_scale(
        self,
        settings,
    ):
        """P2 field T(x): the time scale inside the region, 0 outside, smooth ramps."""
        felics_mesh = self._mean_flow._mesh
        space = fem.functionspace(felics_mesh.dolfinxMesh, ("Lagrange", 2))
        field = fem.Function(space)
        width = float(settings["RampWidth"])
        time_scale = self._time_scale
        if self._region == "Wake":
            x_edge = float(settings["TrailingEdgeX"])

            def shape(points):
                return time_scale * smoothstep((points[0] - (x_edge - width)) / (2.0 * width))
        else:
            x_start, x_end = (float(value) for value in settings["XRange"])

            def shape(points):
                return time_scale * (
                    smoothstep((points[0] - (x_start - width)) / (2.0 * width))
                    * (1.0 - smoothstep((points[0] - (x_end - width)) / (2.0 * width)))
                )
        field.interpolate(shape)
        values = np.real(field.x.array)
        logger.info(
            f"  lag time scale restricted to region {self._region}: "
            f"{np.count_nonzero(values > 1e-12)} of {values.size} P2 values nonzero."
        )
        return field

    # ------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------
    @property
    def active(self):
        """True if the eddy viscosity responds with a lag."""
        return self._model != "None" and self._time_scale > 0.0

    @property
    def eigenvalue(self):
        """Eigenvalue used in the frequency-domain factor."""
        return self._eigenvalue

    @property
    def model(self):
        """Name of the lag model."""
        return self._model

    @property
    def time_scale(self):
        """Time scale T (uniform value; inside the region if one is set)."""
        return self._time_scale

    def factor(self):
        """Complex factor multiplying the eddy viscosity seen by the perturbation.

        Returns
        -------
        ufl expression or complex
            mu_t / (1 - i lambda T) ("Relaxation") or mu_t exp(i lambda T)
            ("Delay") divided by mu_t; 1 if the model is inactive.
        """
        if not self.active:
            return 1.0
        time_scale = self._time_scale_field if self._time_scale_field is not None else self._time_scale
        argument = 1j * self._eigenvalue * time_scale
        if self._model == "Delay":
            return ufl.exp(argument) if self._time_scale_field is not None else complex(np.exp(argument))
        return 1.0 / (1.0 - argument)
