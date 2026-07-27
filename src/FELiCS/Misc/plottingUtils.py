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

# Third party libraries
import  dolfinx
import  basix.ufl                   as ufl
import  matplotlib.pyplot           as plt
from    matplotlib.tri              import Triangulation
from    mpl_toolkits.axes_grid1     import make_axes_locatable
import  numpy                       as np

# Local libraries
from    FELiCS.Misc.logging         import Logger
logger = Logger.get_logger("felics")


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _plot_domain_boundaries(dolfinx_mesh, axes, tdim):
    """Plot the exterior boundary edges of *dolfinx_mesh* as black lines on *axes*.

    Parameters
    ----------
    dolfinx_mesh : dolfinx.mesh.Mesh
    axes : matplotlib.axes.Axes
    tdim : int
        Topological dimension of the mesh.
    """
    coords              = dolfinx_mesh.geometry.x
    boundary_facets     = dolfinx.mesh.exterior_facet_indices(dolfinx_mesh.topology)
    dolfinx_mesh.topology.create_connectivity(tdim - 1, 0)
    facet_to_vertices   = dolfinx_mesh.topology.connectivity(tdim - 1, 0)
    for facet_idx in boundary_facets:
        start           = facet_to_vertices.offsets[facet_idx]
        end             = facet_to_vertices.offsets[facet_idx + 1]
        vertices        = facet_to_vertices.array[start:end]
        edge_coords     = coords[vertices]
        axes.plot(
            edge_coords[:, 0],
            edge_coords[:, 1],
            'k-',
            linewidth=0.8,
            alpha=0.5,
        )


def _add_colorbar_for_contour(axes, contour, label):
    """Attach a narrow colorbar to *axes* for *contour*.

    Parameters
    ----------
    axes    : matplotlib.axes.Axes
    contour : matplotlib.collections.TriMesh or similar mappable
    label   : str
    """
    divider      = make_axes_locatable(axes)
    colorbar_axes = divider.append_axes("right", size="2%", pad=0.5)
    cbar = plt.colorbar(contour, label=label, cax=colorbar_axes)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()

# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def plot_modal_spectrum(eigval, eigval_adjoint=None, ax=None, xlim=None, ylim=None):
    """Plot direct and optional adjoint eigenvalues.

    Parameters
    ----------
    eigval : numpy.ndarray
        Complex direct eigenvalues.
    eigval_adjoint : numpy.ndarray, optional
        Complex adjoint eigenvalues.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created when *None*.
    xlim : tuple, optional
        Limits for the x-axis.
    ylim : tuple, optional
        Limits for the y-axis.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    stable   = eigval[eigval.imag <= 0]
    unstable = eigval[eigval.imag > 0]

    if len(stable) > 0:
        ax.scatter(stable.real, stable.imag, s=20, alpha=0.6,
                   c='blue', label=f'Stable ({len(stable)})', edgecolors='none')
    if len(unstable) > 0:
        ax.scatter(unstable.real, unstable.imag, s=30, alpha=0.8,
                   c='red', label=f'Unstable ({len(unstable)})',
                   edgecolors='black', linewidths=0.5)
    if eigval_adjoint is not None and len(eigval_adjoint) > 0:
        ax.scatter(eigval_adjoint.real, eigval_adjoint.imag, s=20, alpha=0.6,
                   c='green', label=f'Adjoint ({len(eigval_adjoint)})', marker='x')

    ax.set_xlabel(r'$\mathrm{Re}(\omega)$', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Im}(\omega)$', fontsize=12)
    ax.set_title('Eigenvalue Spectrum (Modal Analysis)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax.legend(fontsize=10)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    plt.tight_layout()
    return ax


def plot_resolvent_spectrum(omega, gains, ax=None, xlim=None, ylim=None):
    """Plot resolvent gain curves.
    
    Parameters
    ----------
    omega : numpy.ndarray
        Frequencies at which gains are computed.
    gains : numpy.ndarray
        Gains corresponding to each frequency. Should be a 2D array where each column represents a different mode.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created when *None*.
    xlim : tuple, optional
        Limits for the x-axis.
    ylim : tuple, optional
        Limits for the y-axis.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.hot(np.linspace(0, 1, gains.shape[1]))
    for i in range(gains.shape[1]):
        ax.plot(omega, gains[:, i], 'o-', label=f'Mode {i}',
                color=colors[i], linewidth=2, markersize=4, alpha=0.7)

    ax.set_xlabel(r'$\omega$', fontsize=12)
    ax.set_ylabel(r'Gains squared $\sigma^2$', fontsize=12)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    plt.tight_layout()
    return ax


# TODO: check responsibilities
# Sophies first idea:
# - give only dolfinx mesh, scalar function, and all the plot info
#   (NOT field, felicsMesh etc)

def plot_field(
    field,
    variableName      = None,
    xlim              = None,
    ylim              = None,
    plotType          = "real",
    clim              = None,
    axes              = None,
    showBoundaries    = True,
    free_aspect_ratio = False,
):
    """Plot a scalar (or scalar component of a vector/mixed) FELiCS field.

    This is a debug-level visualisation utility, not intended for
    publication-quality figures.

    Parameters
    ----------
    field             : FELiCS Field
    variableName      : str, optional
        Subfield or component name to plot (e.g. ``'u_x'``).
    xlim              : tuple, optional
    ylim              : tuple, optional
    plotType          : {'real', 'imag', 'magnitude', 'angle'}
    clim              : tuple, optional
    axes              : matplotlib.axes.Axes, optional
    showBoundaries    : bool
    free_aspect_ratio : bool

    Returns
    -------
    matplotlib.axes.Axes
    """
    if field.mesh.dolfinxMesh.topology.dim != 2:
        logger.warning(
            "Field.plot(): plotting is currently only implemented for 2D "
            "meshes. Returning without plotting."
        )
        return

    # TODO: maybe just use the selected field instead of re-calling the method?
    # Delegate vector / mixed fields to the scalar component
    if field.space.num_sub_spaces > 1:
        selected = field._select_field_for_plot(field, variableName)
        return plot_field(
            selected,
            variableName      = variableName,
            xlim              = xlim,
            ylim              = ylim,
            plotType          = plotType,
            clim              = clim,
            axes              = axes,
            showBoundaries    = showBoundaries,
            free_aspect_ratio = free_aspect_ratio,
        )

    mesh = field.mesh.dolfinxMesh
    u_h  = field.function
    tdim = mesh.topology.dim

    mesh.topology.create_connectivity(tdim, 0)

    # Interpolate to P1 so values sit at vertices
    V1 = dolfinx.fem.functionspace(mesh, ufl.element("CG", "triangle", 1))
    u1 = dolfinx.fem.Function(V1)
    u1.interpolate(u_h)

    if plotType == "imag":
        phi_vertex = np.imag(u1.x.array)
        cmap       = "seismic"
    elif plotType == "magnitude":
        phi_vertex = np.abs(u1.x.array)
        cmap       = "magma"
    elif plotType == "angle":
        phi_vertex = np.angle(u1.x.array)
        cmap       = "hsv"
    else:
        phi_vertex = np.real(u1.x.array)
        cmap       = "seismic"

    cells_to_vertices = mesh.topology.connectivity(tdim, 0).array
    triangles         = cells_to_vertices.reshape(-1, 3)
    coords            = mesh.geometry.x
    triang            = Triangulation(coords[:, 0], coords[:, 1], triangles=triangles)

    if axes is None:
        _, axes = plt.subplots()

    if clim is None:
        if plotType == "magnitude":
            clim = (0, np.max(phi_vertex))
        else:
            clim = (-0.5 * np.max(np.abs(phi_vertex)), 0.5 * np.max(np.abs(phi_vertex)))

    contour = axes.tripcolor(
        triang,
        phi_vertex,
        shading = 'gouraud',
        cmap    = cmap,
        vmin    = clim[0],
        vmax    = clim[1],
    )

    if showBoundaries:
        _plot_domain_boundaries(mesh, axes, tdim)

    axes.set_xlabel('x')
    axes.set_ylabel('y')

    title = field.name if field.name else "scalar_field"
    if plotType == "imag":
        title += "_imag"
    elif plotType == "magnitude":
        title += "_magnitude"
    elif plotType == "angle":
        title += "_angle"
    else:
        title += "_real"

    if not free_aspect_ratio:
        axes.set_aspect('equal')
    axes.grid(True, alpha=0.3)

    _add_colorbar_for_contour(axes, contour, title)

    if xlim is not None:
        axes.set_xlim(xlim)
    if ylim is not None:
        axes.set_ylim(ylim)

    plt.tight_layout()
    return axes
