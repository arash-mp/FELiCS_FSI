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
Unit tests of the rigid-body FSI coupling (FSI.MovingBaseState option and the
per-mode aerodynamic force diagnostics).

All tests are validation tests against references implemented here:

* bordering of A with the moving-base-state columns  vs  a dense construction
* static check                                       vs  a dense solve
* aerodynamic_force_diagnostics                      vs  the 1-DOF formulas
* weak form of the momentum forcing                  vs  the strong form + the
                                                         integration-by-parts
                                                         boundary terms
* force-row matrices G and Geo                       vs  analytical integrals
* selection of the force-row terms (MovingBaseStateForceTerms)
"""
# Standard libraries
import os
import sys
from types import SimpleNamespace

# Third party libraries
import numpy as np
import pytest
from mpi4py import MPI
from petsc4py import PETSc
from basix.ufl import element, mixed_element
from dolfinx import mesh as dolfinx_mesh
from dolfinx.fem import Expression, Function, functionspace, form, assemble_scalar
from dolfinx.fem.petsc import LinearProblem
from ufl import (
    FacetNormal,
    TestFunction,
    TrialFunction,
    conditional,
    gt,
    lt,
    Measure,
    SpatialCoordinate,
    TestFunctions,
    as_vector,
    div,
    dot,
    grad,
    inner,
    sym,
)

# Local libraries
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from FELiCS.Equation.FSI.RigidBodyFSI import RigidBodyMotionFSI, _RigidBodyDOF
from FELiCS.Misc.tensorUtils import CoordinateSystem, Tensor


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _bare_fsi(n_dof, n_fluid, M, D, K):
    """A RigidBodyMotionFSI with only the attributes the algebraic methods
    need (no FELiCS case set-up)."""
    fsi = RigidBodyMotionFSI.__new__(RigidBodyMotionFSI)
    fsi.n_dof = n_dof
    fsi.n_struct = 2 * n_dof
    fsi.N_fluid = n_fluid
    fsi.M_mat = np.asarray(M, dtype=float)
    fsi.D_mat = np.asarray(D, dtype=float)
    fsi.K_mat = np.asarray(K, dtype=float)
    fsi.moving_base_state = True
    fsi._mbs_columns = None
    fsi._mbs_force_matrix = np.zeros((n_dof, n_dof), dtype=complex)
    fsi._force_arrays = None
    return fsi


def _petsc_from_dense(dense):
    rows, cols = np.nonzero(dense)
    n = dense.shape[0]
    mat = PETSc.Mat().createAIJ((n, n), nnz=(np.bincount(rows, minlength=n) + 1).astype(PETSc.IntType),
                                comm=PETSc.COMM_SELF)
    for r, c in zip(rows, cols):
        mat.setValue(int(r), int(c), dense[r, c])
    mat.assemble()
    return mat


def _random_case(seed=0, n_fluid=40, n_dof=2):
    rng = np.random.default_rng(seed)
    dense = (rng.standard_normal((n_fluid, n_fluid))
             + 1j * rng.standard_normal((n_fluid, n_fluid)))
    dense[rng.random((n_fluid, n_fluid)) < 0.7] = 0.0
    dense += 10.0 * np.eye(n_fluid)
    interface = np.array([3, 7, 8, 21])
    couplings = [
        (int(row),
         rng.standard_normal(n_dof) + 1j * rng.standard_normal(n_dof),
         rng.standard_normal(n_dof))
        for row in interface
    ]
    force = [rng.standard_normal(n_fluid) * (rng.random(n_fluid) < 0.5)
             + 0j for _ in range(n_dof)]
    columns = [rng.standard_normal(n_fluid) * (rng.random(n_fluid) < 0.6)
               + 0j for _ in range(n_dof)]
    C = rng.standard_normal((n_dof, n_dof)) + 0j
    M = np.array([[2.0, 0.3], [0.3, 1.5]])[:n_dof, :n_dof]
    Dm = np.diag([0.1, 0.2])[:n_dof, :n_dof]
    Km = np.diag([3.0, 4.0])[:n_dof, :n_dof]
    return dense, interface, couplings, force, columns, C, M, Dm, Km


def _dense_reference_operator(dense, interface, couplings, force, columns, C, K, D):
    """Augmented A written down directly from the documented layout."""
    n_fluid = dense.shape[0]
    n_dof = len(force)
    A = np.zeros((n_fluid + 2 * n_dof, n_fluid + 2 * n_dof), dtype=complex)
    A[:n_fluid, :n_fluid] = dense
    eta = lambda k: n_fluid + 2 * k
    phi = lambda k: n_fluid + 2 * k + 1
    for j in range(n_dof):
        mask = columns[j] != 0
        mask[interface] = False
        A[:n_fluid, eta(j)] += 1j * columns[j] * mask
    for row, b_list, a_list in couplings:
        A[row, :] = 0.0
        A[row, row] = 1.0
        for k in range(n_dof):
            A[row, eta(k)] = b_list[k]
            A[row, phi(k)] = -a_list[k]
    for k in range(n_dof):
        A[eta(k), phi(k)] = 1j
        A[phi(k), :n_fluid] = 1j * force[k]
        for j in range(n_dof):
            A[phi(k), eta(j)] = -1j * K[k, j] + 1j * C[k, j]
            A[phi(k), phi(j)] = -1j * D[k, j]
    return A


# ---------------------------------------------------------------------------
# algebra: bordering, static check, diagnostics
# ---------------------------------------------------------------------------
def test_border_operator_with_moving_base_state_columns():
    dense, interface, couplings, force, columns, C, M, Dm, Km = _random_case()
    fsi = _bare_fsi(2, dense.shape[0], M, Dm, Km)
    fsi._mbs_columns = columns
    fsi._mbs_force_matrix = C
    A_ff = _petsc_from_dense(dense)
    force_vecs = [PETSc.Vec().createWithArray(f.copy(), comm=PETSc.COMM_SELF)
                  for f in force]
    A_aug = fsi._border_A(A_ff, force_vecs, couplings)
    ai, aj, av = A_aug.getValuesCSR()
    n = A_aug.getSize()[0]
    result = np.zeros((n, n), dtype=complex)
    for r in range(n):
        result[r, aj[ai[r]:ai[r + 1]]] = av[ai[r]:ai[r + 1]]
    reference = _dense_reference_operator(dense, interface, couplings, force,
                                   columns, C, Km, Dm)
    assert np.allclose(result, reference, rtol=0.0, atol=1e-14)
    # CSR column indices must be sorted in every row
    for r in range(n):
        assert np.all(np.diff(aj[ai[r]:ai[r + 1]]) > 0)


def test_static_check_matches_dense_solve():
    dense, interface, couplings, force, columns, C, M, Dm, Km = _random_case(1)
    n_fluid = dense.shape[0]
    fsi = _bare_fsi(2, n_fluid, M, Dm, Km)
    fsi._mbs_columns = columns
    fsi._mbs_force_matrix = C
    fsi._force_arrays = force
    A_ff = _petsc_from_dense(dense)
    force_vecs = [PETSc.Vec().createWithArray(f.copy(), comm=PETSc.COMM_SELF)
                  for f in force]
    A_aug = fsi._border_A(A_ff, force_vecs, couplings)
    reference = _dense_reference_operator(dense, interface, couplings, force,
                                   columns, C, Km, Dm)
    expected = np.zeros((2, 2), dtype=complex)
    for j in range(2):
        q = np.linalg.solve(reference[:n_fluid, :n_fluid],
                            -reference[:n_fluid, n_fluid + 2 * j])
        for k in range(2):
            expected[k, j] = force[k] @ q + C[k, j]

    result = fsi._static_force_matrices(A_aug, couplings)
    assert np.allclose(result["moving base state"], expected,
                       rtol=1e-10, atol=1e-12)
    frozen = np.zeros((2, 2), dtype=complex)
    for j in range(2):
        rhs = np.zeros(n_fluid, dtype=complex)
        for row, b_list, a_list in couplings:
            rhs[row] = -b_list[j]
        q = np.linalg.solve(reference[:n_fluid, :n_fluid], rhs)
        for k in range(2):
            frozen[k, j] = force[k] @ q
    assert np.allclose(result["frozen base state (original formulation)"],
                       frozen, rtol=1e-10, atol=1e-12)


def test_aerodynamic_force_diagnostics_single_dof():
    mass, damping, stiffness = 2.5, 0.3, 4.0
    fsi = _bare_fsi(1, 3, [[mass]], [[damping]], [[stiffness]])
    fsi._force_arrays = [np.array([1.0, 0.0, 0.0], dtype=complex)]
    fsi._mbs_force_matrix = np.array([[0.7]], dtype=complex)
    eigen_value = 0.83 + 0.041j
    s = -1j * eigen_value
    eta = 0.4 - 0.9j
    phi = s * eta
    total_force = mass * s * phi + damping * phi + stiffness * eta
    # the fluid part carries the force minus the moving-base-state part
    q_aug = np.array([total_force - 0.7 * eta, 0.0, 0.0, eta, phi])
    diag = fsi.aerodynamic_force_diagnostics(eigen_value, q_aug)
    omega_r, omega_i = eigen_value.real, eigen_value.imag
    assert diag["residual"] < 1e-13
    assert diag["phi_minus_s_eta"] < 1e-13
    assert np.isclose(diag["k_aero"][0],
                      mass * (omega_r**2 + omega_i**2) - stiffness)
    assert np.isclose(diag["c_aero"][0], -2.0 * mass * omega_i - damping)


# ---------------------------------------------------------------------------
# finite-element forms
# ---------------------------------------------------------------------------
def _fe_case():
    domain = dolfinx_mesh.create_rectangle(
        MPI.COMM_WORLD,
        ((0.0, 0.0), (1.0, 0.8)),
        (6, 5),
        cell_type=dolfinx_mesh.CellType.triangle,
    )
    cell = domain.topology.cell_type.name
    p2_vector = element("Lagrange", cell, 2, shape=(2,))
    p1 = element("Lagrange", cell, 1)
    mixed = functionspace(domain, mixed_element([p2_vector, p1]))
    vector_space = functionspace(domain, p2_vector)
    scalar_space = functionspace(domain, ("Lagrange", 2))

    fdim = domain.topology.dim - 1
    domain.topology.create_connectivity(fdim, domain.topology.dim)
    facets = dolfinx_mesh.exterior_facet_indices(domain.topology)
    tags = dolfinx_mesh.meshtags(domain, fdim, facets,
                                 np.full(len(facets), 1, dtype=np.int32))
    x = SpatialCoordinate(domain)
    cs = CoordinateSystem(x, "cartesian", 2)
    test_u, test_p = TestFunctions(mixed)
    return SimpleNamespace(
        domain=domain, mixed=mixed, vector_space=vector_space,
        scalar_space=scalar_space, x=x, cs=cs,
        ds_tagged=Measure("ds", domain=domain, subdomain_data=tags),
        ds_all=Measure("ds", domain=domain),
        test=[Tensor(test_u, cs, mayHaveSpectralDimension=True),
              Tensor(test_p, cs, mayHaveSpectralDimension=True)],
    )


def _fe_fsi(case, mean_u, mean_mu):
    fsi = RigidBodyMotionFSI.__new__(RigidBodyMotionFSI)
    fsi.cs = case.cs
    fsi.x = case.x
    fsi.n = Tensor(FacetNormal(case.domain), case.cs)
    fsi.X = case.test
    fsi.u_index = 0
    fsi.mesh = case.domain
    fsi.interface_id = 1
    fsi.moving_base_state = True
    fsi.mbs_wake = False
    fsi.mbs_wake_length = 0.5
    fsi.mbs_wake_rotation = "Translate"
    fsi.mbs_wake_translation = "Rigid"
    fsi._trailing_edge = None
    fsi.FEMSpaces = SimpleNamespace(VMixed=case.mixed)
    fsi.meanFlow = SimpleNamespace(
        u=Tensor(mean_u, case.cs),
        mu_tot=Tensor(mean_mu, case.cs),
    )
    fsi.equation = SimpleNamespace(
        all_ds=case.ds_all,
        ds=case.ds_tagged,
        boundaryHandler=SimpleNamespace(
            get_list_of_dirichlet_b_cs_for_dolfinx=lambda space, exclude_ids: [],
        ),
    )
    return fsi


def _interpolate(space, callable_):
    function = Function(space)
    function.interpolate(callable_)
    return function


def test_moving_base_state_vector_matches_strong_form():
    case = _fe_case()
    # quadratic U0 and mu are represented exactly by P2, so S0 and grad(mu)
    # are continuous and the integration by parts has no interior-facet
    # jump terms (the residual may be any P2 field)
    mean_u = _interpolate(case.vector_space, lambda x: np.vstack((
        0.2 + x[0] * x[1] - 0.3 * x[1]**2, 0.5 * x[0]**2 - x[0] * x[1] + 0.1)))
    mean_mu = _interpolate(case.scalar_space,
                           lambda x: 0.01 * (1.0 + x[0]**2 + 0.5 * x[0] * x[1]))
    residual = _interpolate(case.vector_space, lambda x: np.vstack((
        x[0] * x[1] - 0.3 * x[1]**2, np.exp(0.4 * x[0]) - x[1])))
    # w: a velocity field that does NOT vanish on the boundary, p-part zero
    w_mixed = Function(case.mixed)
    w_mixed.sub(0).interpolate(lambda x: np.vstack((
        1.0 + x[0] * x[1], x[0] - 0.7 * x[1]**2)))
    w = w_mixed.sub(0)
    fsi = _fe_fsi(case, mean_u, mean_mu)
    n = FacetNormal(case.domain)
    S0 = sym(grad(mean_u))

    dofs = [
        _RigidBodyDOF("translation", 1.0, 0.0, 0.0, direction=(0.0, 1.0)),
        _RigidBodyDOF("translation", 1.0, 0.0, 0.0, direction=(0.6, 0.8)),
        _RigidBodyDOF("rotation", 1.0, 0.0, 0.0, center=(0.41, 0.1)),
    ]
    for dof in dofs:
        d = dof.displacement_ufl(case.x)
        rotation = 1.0 if dof.motion_type == "rotation" else 0.0
        g = fsi._moving_base_state_vector(dof, Tensor(residual, case.cs))
        assembled = np.dot(g, w_mixed.x.array)

        # strong form + the boundary terms of the integration by parts that
        # the weak form keeps (convective) / drops (viscous)
        m = dot(grad(mean_mu), d)
        z_cross_r = as_vector((-residual[1], residual[0]))
        f = -grad(residual) * d + rotation * z_cross_r - div(2 * m * S0)
        reference = assemble_scalar(form(
            dot(f, w) * Measure("dx", domain=case.domain)
            + dot(2 * m * S0 * n, w) * case.ds_all
        ))
        assert np.isclose(assembled, reference, rtol=1e-10, atol=1e-13), \
            (dof.label, assembled, reference)


def test_moving_base_state_force_matrices_analytical():
    case = _fe_case()
    # U0 = (y, x): S0 = [[0, 1], [1, 0]];  mu = x^2/2: grad mu = (x, 0)
    mean_u = _interpolate(case.vector_space,
                          lambda x: np.vstack((x[1], x[0])))
    mean_mu = _interpolate(case.scalar_space, lambda x: 0.5 * x[0]**2)
    residual = _interpolate(case.vector_space,
                            lambda x: np.vstack((x[0], x[1])))
    fsi = _fe_fsi(case, mean_u, mean_mu)
    fsi.dofs = [
        _RigidBodyDOF("translation", 1.0, 0.0, 0.0, direction=(1.0, 0.0)),
        _RigidBodyDOF("translation", 1.0, 0.0, 0.0, direction=(0.0, 1.0)),
    ]
    fsi.n_dof = 2
    G, Geo = fsi._moving_base_state_force_matrices(Tensor(residual, case.cs))
    # closed contour, n = outward normal of the domain, area 0.8:
    #   Geo_kj = oint (d_j.n)(d_k.R0) = int div((d_k.R0) d_j) = 0.8 delta_kj
    #   G_kj   = oint 2 (d_j.grad mu) d_k.(S0 n) = int div(2 (d_j.grad mu) S0 d_k)
    #          = 2 (S0 d_k)_x * 0.8 for j = x, 0 for j = y
    area = 0.8
    assert np.allclose(Geo, area * np.eye(2), atol=1e-12)
    assert np.allclose(G, np.array([[0.0, 0.0], [2.0 * area, 0.0]]), atol=1e-12)


def test_selected_mbs_force_matrix():
    fsi = _bare_fsi(2, 4, np.eye(2), np.zeros((2, 2)), np.eye(2))
    G = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=complex)
    Geo = np.array([[10.0, 20.0], [30.0, 40.0]], dtype=complex)
    expected = {"All": G + Geo, "NoGeo": G, "NoG": Geo, "None": np.zeros((2, 2))}
    for terms, reference in expected.items():
        fsi.mbs_force_terms = terms
        assert np.allclose(fsi._selected_mbs_force_matrix(G, Geo), reference), terms
    # default (attribute not set, e.g. older pickles): all terms
    del fsi.mbs_force_terms
    assert np.allclose(fsi._selected_mbs_force_matrix(G, Geo), G + Geo)


# ---------------------------------------------------------------------------
# rotation with the wake-following transport (FSI.MovingBaseStateWake)
# ---------------------------------------------------------------------------
_WAKE_START = 1.0 / 3.0      # trailing edge x; with the length below the blend
_WAKE_LENGTH = 0.5           # ends at x = 5/6: both on mesh lines (6 cells)


def _reference_wake_transport(x, center, trailing_edge, length, frozen=False):
    """Independent implementation of the transport field (conditional
    instead of max/min, explicit smoothstep).  frozen: the wake is not moved
    (WakeRotation = "Frozen") instead of translated by the TE displacement."""
    x_te, y_te = trailing_edge
    cx, cy = center
    ramp = (x[0] - x_te) / length
    ramp = conditional(lt(ramp, 0.0), 0.0, conditional(gt(ramp, 1.0), 1.0, ramp))
    weight = 1.0 - ramp * ramp * (3.0 - 2.0 * ramp)
    rigid = as_vector((-(x[1] - cy), x[0] - cx))
    if frozen:
        return weight * rigid
    shifted = as_vector((-(y_te - cy), x_te - cx))
    return weight * rigid + (1.0 - weight) * shifted


@pytest.mark.parametrize("wake_rotation", ["Translate", "Frozen"])
def test_wake_transport_vector_matches_strong_form(wake_rotation):
    case = _fe_case()
    mean_u = _interpolate(case.vector_space, lambda x: np.vstack((
        0.2 + x[0] * x[1] - 0.3 * x[1]**2, 0.5 * x[0]**2 - x[0] * x[1] + 0.1)))
    mean_mu = _interpolate(case.scalar_space,
                           lambda x: 0.01 * (1.0 + x[0]**2 + 0.5 * x[0] * x[1]))
    residual = _interpolate(case.vector_space, lambda x: np.vstack((
        x[0] * x[1] - 0.3 * x[1]**2, np.exp(0.4 * x[0]) - x[1])))
    w_mixed = Function(case.mixed)
    w_mixed.sub(0).interpolate(lambda x: np.vstack((
        1.0 + x[0] * x[1], x[0] - 0.7 * x[1]**2)))
    w = w_mixed.sub(0)
    fsi = _fe_fsi(case, mean_u, mean_mu)
    fsi.mbs_wake = True
    fsi.mbs_wake_length = _WAKE_LENGTH
    fsi.mbs_wake_rotation = wake_rotation
    fsi._trailing_edge = (_WAKE_START, 0.1)
    dof = _RigidBodyDOF("rotation", 1.0, 0.0, 0.0, center=(0.2, 0.15))

    g = fsi._moving_base_state_vector(dof, Tensor(residual, case.cs))
    assembled = np.dot(g, w_mixed.x.array)

    # strong form f = -(xi.grad)R0 + a z x R0 - div(2 (xi.grad mu) S0),
    # a = (curl xi)/2, plus the dropped viscous boundary term
    xi = _reference_wake_transport(case.x, dof.center, fsi._trailing_edge,
                                   _WAKE_LENGTH, frozen=(wake_rotation == "Frozen"))
    rotation = 0.5 * (xi[1].dx(0) - xi[0].dx(1))
    n = FacetNormal(case.domain)
    S0 = sym(grad(mean_u))
    m = dot(grad(mean_mu), xi)
    z_cross_r = as_vector((-residual[1], residual[0]))
    f = -grad(residual) * xi + rotation * z_cross_r - div(2 * m * S0)
    reference = assemble_scalar(form(
        dot(f, w) * Measure("dx", domain=case.domain)
        + dot(2 * m * S0 * n, w) * case.ds_all
    ))
    assert np.isclose(assembled, reference, rtol=1e-10, atol=1e-13), \
        (assembled, reference)


def test_wake_transport_limits():
    case = _fe_case()
    fsi = _fe_fsi(case, Function(case.vector_space), Function(case.scalar_space))
    fsi.mbs_wake = True
    fsi.mbs_wake_length = _WAKE_LENGTH
    fsi._trailing_edge = (_WAKE_START, 0.1)
    center = (0.2, 0.15)
    rotation_dof = _RigidBodyDOF("rotation", 1.0, 0.0, 0.0, center=center)
    xi, rotation, divergence = fsi._transport_field(rotation_dof)
    assert rotation is not None and divergence is not None

    xi_h = Function(case.vector_space)
    xi_h.interpolate(Expression(
        xi, case.vector_space.element.interpolation_points()))
    coords = case.vector_space.tabulate_dof_coordinates()
    values = xi_h.x.array.real.reshape(-1, 2)
    rigid = np.c_[-(coords[:, 1] - center[1]), coords[:, 0] - center[0]]
    shifted = np.array([-(0.1 - center[1]), _WAKE_START - center[0]])
    upstream = coords[:, 0] <= _WAKE_START + 1e-12
    downstream = coords[:, 0] >= _WAKE_START + _WAKE_LENGTH - 1e-12
    assert upstream.any() and downstream.any()
    assert np.allclose(values[upstream], rigid[upstream], atol=1e-12)
    assert np.allclose(values[downstream], shifted, atol=1e-12)

    # translations are not affected by the option
    translation_dof = _RigidBodyDOF("translation", 1.0, 0.0, 0.0,
                                    direction=(0.0, 1.0))
    xi_t, rotation_t, divergence_t = fsi._transport_field(translation_dof)
    assert rotation_t is None and divergence_t is None

    # without the option a rotation is rigid everywhere
    fsi.mbs_wake = False
    xi_r, rotation_r, divergence_r = fsi._transport_field(rotation_dof)
    assert rotation_r == 1.0 and divergence_r is None


def _interpolated_transport(case, xi):
    xi_h = Function(case.vector_space)
    xi_h.interpolate(Expression(
        xi, case.vector_space.element.interpolation_points()))
    return xi_h.x.array.real.reshape(-1, 2)


def test_wake_transport_frozen_and_rigid_modes():
    """WakeRotation "Frozen" / "Rigid" and WakeTranslation "Frozen" / "Rigid":
    rigid upstream of the TE, the chosen wake transport downstream."""
    case = _fe_case()
    fsi = _fe_fsi(case, Function(case.vector_space), Function(case.scalar_space))
    fsi.mbs_wake = True
    fsi.mbs_wake_length = _WAKE_LENGTH
    fsi._trailing_edge = (_WAKE_START, 0.1)
    center = (0.2, 0.15)
    coords = case.vector_space.tabulate_dof_coordinates()
    upstream = coords[:, 0] <= _WAKE_START + 1e-12
    downstream = coords[:, 0] >= _WAKE_START + _WAKE_LENGTH - 1e-12
    assert upstream.any() and downstream.any()
    rotation_dof = _RigidBodyDOF("rotation", 1.0, 0.0, 0.0, center=center)
    heave_dof = _RigidBodyDOF("translation", 1.0, 0.0, 0.0, direction=(0.0, 1.0))
    rigid_rot = np.c_[-(coords[:, 1] - center[1]), coords[:, 0] - center[0]]
    heave = np.array([0.0, 1.0])

    # pitch: wake not moved
    fsi.mbs_wake_rotation = "Frozen"
    xi, rotation, divergence = fsi._transport_field(rotation_dof)
    assert rotation is not None and divergence is not None
    values = _interpolated_transport(case, xi)
    assert np.allclose(values[upstream], rigid_rot[upstream], atol=1e-12)
    assert np.allclose(values[downstream], 0.0, atol=1e-12)

    # pitch: rigid everywhere, as without the option
    fsi.mbs_wake_rotation = "Rigid"
    xi, rotation, divergence = fsi._transport_field(rotation_dof)
    assert rotation == 1.0 and divergence is None

    # heave: follows the body everywhere by default
    xi, rotation, divergence = fsi._transport_field(heave_dof)
    assert rotation is None and divergence is None

    # heave: wake not moved
    fsi.mbs_wake_translation = "Frozen"
    xi, rotation, divergence = fsi._transport_field(heave_dof)
    assert rotation is not None and divergence is not None
    values = _interpolated_transport(case, xi)
    assert np.allclose(values[upstream], heave, atol=1e-12)
    assert np.allclose(values[downstream], 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# discrete base-flow residual (FSI.MovingBaseStateResidual = "Discrete")
# ---------------------------------------------------------------------------
def test_discrete_base_residual_matches_projected_strong_form():
    case = _fe_case()
    # quadratic U0, P0 and mu: grad U0 and the stress are continuous, so the
    # discrete residual has no inter-element jump terms and must equal the L2
    # projection of the strong residual (grad U0) U0 + grad P0 - div tau0
    mean_u = _interpolate(case.vector_space, lambda x: np.vstack((
        0.2 + x[0] * x[1] - 0.3 * x[1]**2, 0.5 * x[0]**2 - x[0] * x[1] + 0.1)))
    mean_p = _interpolate(case.scalar_space,
                          lambda x: 0.3 * x[0]**2 - 0.2 * x[0] * x[1] + x[1])
    mean_mu = _interpolate(case.scalar_space,
                           lambda x: 0.01 * (1.0 + x[0]**2 + 0.5 * x[0] * x[1]))
    density = _interpolate(case.scalar_space,
                           lambda x: np.ones(x.shape[1]))
    fsi = _fe_fsi(case, mean_u, mean_mu)
    fsi.meanFlow.p = Tensor(mean_p, case.cs)
    fsi.meanFlow.rho = Tensor(density, case.cs)
    residual = fsi._base_residual_discrete().ufl_tens

    tau = mean_mu * (grad(mean_u) + grad(mean_u).T)
    strong = grad(mean_u) * mean_u + grad(mean_p) - div(tau)
    trial = TrialFunction(case.vector_space)
    test = TestFunction(case.vector_space)
    reference = LinearProblem(
        inner(trial, test) * Measure("dx", domain=case.domain),
        inner(strong, test) * Measure("dx", domain=case.domain),
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    ).solve()
    assert np.allclose(residual.x.array, reference.x.array, rtol=1e-9, atol=1e-11)
    assert np.abs(reference.x.array).max() > 1e-2
