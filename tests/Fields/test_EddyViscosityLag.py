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
Unit tests of the lagged eddy viscosity (Case.EddyViscosityLag).  Validation
tests against references implemented here:

* the default and the settings checks (no silent defaults for model keys)
* the frequency-domain factors of the relaxation and delay models
* a zero time scale recovers the frozen model exactly
* the region masks (Wake, Band): time scale inside, zero outside
"""
# Standard libraries
import os
import sys
from types import SimpleNamespace

# Third party libraries
import numpy as np
import pytest
import ufl
from mpi4py import MPI
from dolfinx import fem
from dolfinx import mesh as dolfinx_mesh

# Local libraries
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from FELiCS.Fields.EddyViscosityLag import EddyViscosityLag, smoothstep


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _mean_flow():
    domain = dolfinx_mesh.create_rectangle(
        MPI.COMM_WORLD,
        ((0.0, -0.5), (2.0, 0.5)),
        (40, 10),
        cell_type=dolfinx_mesh.CellType.triangle,
    )
    return SimpleNamespace(_mesh=SimpleNamespace(dolfinxMesh=domain), domain=domain)


def _param(settings):
    return SimpleNamespace(Case=SimpleNamespace(
        EddyViscosityLag=settings,
        CoordinateSystem="Cartesian",
        nDim=2,
    ))


def _uniform(model, time_scale, eigenvalue=(0.7, -0.03)):
    return {"Model": model, "TimeScale": time_scale, "Eigenvalue": list(eigenvalue), "Region": "All"}


def _field_values(lag, mean):
    """Evaluate the (masked) factor at the P2 points, with their x coordinates."""
    space = fem.functionspace(mean.domain, ("Lagrange", 2))
    function = fem.Function(space)
    function.interpolate(fem.Expression(lag.factor(), space.element.interpolation_points()))
    points = space.tabulate_dof_coordinates()
    return points[:, 0], function.x.array


# ---------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------
def test_none_is_default_and_inactive():
    lag = EddyViscosityLag(_param({"Model": "None"}), _mean_flow())
    assert not lag.active
    assert lag.factor() == 1.0


def test_missing_setting_defaults_to_none():
    param = SimpleNamespace(Case=SimpleNamespace(CoordinateSystem="Cartesian", nDim=2))
    assert not EddyViscosityLag(param, _mean_flow()).active


@pytest.mark.parametrize("settings", [
    {"Model": "None", "TimeScale": 1.0},
    {"Model": "Relaxation", "TimeScale": 1.0},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0]},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": 0.7, "Region": "All"},
    {"Model": "Relaxation", "TimeScale": -1.0, "Eigenvalue": [0.7, 0.0], "Region": "All"},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0], "Region": "Body"},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0], "Region": "Band",
     "RampWidth": 0.02},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0], "Region": "Band",
     "RampWidth": 0.02, "XRange": [0.9, 0.5]},
    {"Model": "Relaxation", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0], "Region": "All",
     "RampWidth": 0.02},
    {"Model": "Other", "TimeScale": 1.0, "Eigenvalue": [0.7, 0.0], "Region": "All"},
])
def test_incomplete_or_unknown_settings_are_rejected(settings):
    with pytest.raises(ValueError):
        EddyViscosityLag(_param(settings), _mean_flow())


def test_three_dimensional_case_is_rejected():
    param = SimpleNamespace(Case=SimpleNamespace(
        EddyViscosityLag=_uniform("Relaxation", 1.0), CoordinateSystem="Cartesian", nDim=3))
    with pytest.raises(NotImplementedError):
        EddyViscosityLag(param, _mean_flow())


# ---------------------------------------------------------------------------
# factors
# ---------------------------------------------------------------------------
def test_relaxation_factor_matches_formula():
    eigenvalue = 0.7 - 0.03j
    lag = EddyViscosityLag(_param(_uniform("Relaxation", 2.0, (0.7, -0.03))), _mean_flow())
    assert lag.active
    assert lag.factor() == pytest.approx(1.0 / (1.0 - 1j * eigenvalue * 2.0))


def test_relaxation_lag_and_attenuation():
    """The stress lags the strain by arctan(omega T) and is attenuated by cos of it."""
    omega, time_scale = 0.581, 2.2
    lag = EddyViscosityLag(_param(_uniform("Relaxation", time_scale, (omega, 0.0))), _mean_flow())
    factor = lag.factor()
    angle = np.arctan(omega * time_scale)
    assert np.angle(factor) == pytest.approx(angle)
    assert abs(factor) == pytest.approx(np.cos(angle))


def test_delay_factor_matches_formula():
    eigenvalue = 0.9 + 0.02j
    lag = EddyViscosityLag(_param(_uniform("Delay", 1.7, (0.9, 0.02))), _mean_flow())
    assert lag.factor() == pytest.approx(np.exp(1j * eigenvalue * 1.7))
    assert abs(lag.factor()) == pytest.approx(np.exp(-0.02 * 1.7))


@pytest.mark.parametrize("model", ["Relaxation", "Delay"])
def test_zero_time_scale_is_frozen(model):
    lag = EddyViscosityLag(_param(_uniform(model, 0.0)), _mean_flow())
    assert not lag.active
    assert lag.factor() == 1.0


# ---------------------------------------------------------------------------
# regions
# ---------------------------------------------------------------------------
def test_smoothstep_limits():
    assert smoothstep(np.array([-1.0, 0.0, 0.5, 1.0, 2.0])) == pytest.approx([0.0, 0.0, 0.5, 1.0, 1.0])


def test_band_region_is_time_scale_inside_and_frozen_outside():
    mean = _mean_flow()
    settings = {"Model": "Relaxation", "TimeScale": 2.0, "Eigenvalue": [0.7, 0.0],
                "Region": "Band", "XRange": [0.8, 1.4], "RampWidth": 0.05}
    lag = EddyViscosityLag(_param(settings), mean)
    x, values = _field_values(lag, mean)
    inside = (x > 0.9) & (x < 1.3)
    outside = (x < 0.7) | (x > 1.5)
    assert np.allclose(values[inside], 1.0 / (1.0 - 1j * 0.7 * 2.0))
    assert np.allclose(values[outside], 1.0)


def test_wake_region_is_frozen_upstream_of_the_edge():
    mean = _mean_flow()
    settings = {"Model": "Delay", "TimeScale": 1.5, "Eigenvalue": [0.7, 0.0],
                "Region": "Wake", "TrailingEdgeX": 1.0, "RampWidth": 0.05}
    lag = EddyViscosityLag(_param(settings), mean)
    x, values = _field_values(lag, mean)
    assert np.allclose(values[x < 0.9], 1.0)
    assert np.allclose(values[x > 1.1], np.exp(1j * 0.7 * 1.5))
