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
"""Unit tests for the Numerics.SpongeOnPressure option of SpongeTerm."""

# Standard libraries
import os
import sys
from types import SimpleNamespace

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Third party libraries
import pytest

# Local libraries
import FELiCS.Equation.Equations.SpongeTerm as sponge_module
from FELiCS.Equation.Equations.SpongeTerm import SpongeTerm


class _Absorbing:
    """Stand-in for UFL objects: every product returns the object itself."""

    def __init__(self):
        self.ufl_tens = self

    def __mul__(self, other):
        return self

    __rmul__ = __mul__


class _RecordingForm:
    """Collects the terms passed to ``add``."""

    def __init__(self):
        self.terms = []

    def add(self, term):
        self.terms.append(term)


def _make_sponge(sponge_on_pressure):
    """A SpongeTerm with only the attributes the sponge loop needs."""
    set_of_equations = {
        "Momentum": {"Equation": "NSPrimitive", "Variable": "u"},
        "Mass": {"Equation": "Continuity", "Variable": "p"},
        "Energy": {"Equation": "None", "Variable": "None"},
    }
    term = SpongeTerm.__new__(SpongeTerm)
    term.param = SimpleNamespace(
        Case=SimpleNamespace(
            SetOfEquations=set_of_equations,
            SolutionList=["u", "p"],
        ),
        Numerics=SimpleNamespace(SpongeOnPressure=sponge_on_pressure),
    )
    term.fluc = SimpleNamespace(u=_Absorbing(), p=_Absorbing())
    term.X = [_Absorbing(), _Absorbing()]
    term.J_hat = 1.0
    return term


@pytest.fixture(autouse=True)
def _absorbing_tensor_helpers(monkeypatch):
    monkeypatch.setattr(sponge_module, "i_dot", lambda a, b: _Absorbing())
    monkeypatch.setattr(sponge_module, "i_conj", lambda a: _Absorbing())
    monkeypatch.setattr(sponge_module, "dx", _Absorbing())


@pytest.mark.parametrize(
    "sponge_on_pressure, expected_skipped",
    [(True, False), (False, True)],
)
def test_is_skipped_continuity(sponge_on_pressure, expected_skipped):
    term = _make_sponge(sponge_on_pressure)
    assert term._is_skipped("Mass") is expected_skipped
    assert term._is_skipped("Momentum") is False


@pytest.mark.parametrize(
    "sponge_on_pressure, expected_number_of_terms",
    [(True, 2), (False, 1)],
)
def test_linear_expression_term_count(sponge_on_pressure, expected_number_of_terms):
    term = _make_sponge(sponge_on_pressure)
    weak_form = _RecordingForm()
    mean = SimpleNamespace(spg=_Absorbing())
    term.add_linear_expression(weak_form, mean)
    assert len(weak_form.terms) == expected_number_of_terms
