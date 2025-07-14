import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import matplotlib.tri as tri

# --- User parameters ---
case_path = '/mnt/c/Users/avillie/OneDrive/Documents/PhD/PINNs/c14toWsl/resolvent_tutorial'

gains_file = '/gains.csv'   
mesh_file = '/Resolvent_mesh.h5'        
forcing_file = '/Resolvent_Omega3.1_Forcing_gain0.h5'      
response_file = '/Resolvent_Omega3.1_Response_gain0.h5'      

# --- Load gains and frequencies from CSV ---
df = pd.read_csv(case_path + gains_file)
St = df['omega'] / (2 * np.pi)             # Strouhal number
gains = df.drop(columns=['omega']).values  # shape: (n_freqs, n_resp_modes)
Nmodes = gains.shape[1]                    # Number of solutions computed
# --- Plot gains ---
plt.figure(figsize=(6, 4))
for i in range(Nmodes):
    plt.loglog(St, gains[:, i], label=f'Mode {i+1}', color=plt.cm.Blues(np.linspace(1, 0.3, Nmodes))[i])
plt.xlabel('St')
plt.ylabel('$\\sigma^2$')
plt.title('Resolvent Gains')
plt.legend()
plt.show()

# --- Load and plot a response mode ---

def mesh_triang(mesh):
    mesh_file_handle = h5py.File(mesh, 'r')
    x = np.array(mesh_file_handle['coordinates/x'])
    y = np.array(mesh_file_handle['coordinates/y'])
    connectivity =  mesh_file_handle['cells']['triangles']
    return tri.Triangulation(x, y, triangles=connectivity)
triang = mesh_triang(case_path + mesh_file)

def plot_field(x, y, field, triang, **kwargs):
    '''
    Plot the resolvent modes
    
    Possible extra inputs:
        .ax:            axes handle from existing figure
        .xlims, ylims:  plot limits
        .cmap:          colormap
        .shading:       shading option for tripcolor
    '''
    if not 'ax' in kwargs:
        ax = plt.subplot()
    else:
        ax = kwargs['ax']
        subplot_spec = ax.get_subplotspec()
        total_rows = subplot_spec.get_gridspec().nrows
        current_row = subplot_spec.rowspan.start
        axis3 = not (total_rows is not None and current_row < total_rows - 1)

    if not 'cmap' in kwargs:
        cmap = plt.get_cmap('coolwarm') #'afmhot_r', 'plasma', 'viridis'
    else:
        cmap = kwargs['cmap']
        
    if 'xlims' in kwargs:
        xlims = kwargs['xlims']
        plt.xlim(xlims)
    else: 
        plt.xlim([np.min(x),np.max(x)])
    
    if 'ylims' in kwargs:
        ylims = kwargs['ylims']
        plt.ylim(ylims)
    elif 'triang_r' in kwargs:
        plt.ylim([-np.max(y),np.max(y)])    
    else: 
        plt.ylim([np.min(y),np.max(y)])

    # if 'title' in kwargs:
    #     title = kwargs['title']
        #plt.title(title, x = 0.45, y=0.8)
    title = kwargs['title'] if 'title' in kwargs else ''
    
    if 'shading' in kwargs:
        shading = kwargs['shading']
    else: 
        shading = 'gouraud'#'flat'
    #print(np.min(field), np.max(field), "minmax")
    #vmax, vmin = None,None 
    if ('vmax' in kwargs):
        vmax =  kwargs['vmax']
    else: 
        vmax = None 
    if ('vmin' in kwargs):
        vmin = kwargs['vmin']    
    else: 
        vmin = None    
    if 'triang_r' in kwargs:
        tpc = ax.tripcolor(triang, field, cmap=cmap, vmax = vmax, vmin = vmin, shading=shading)
        triang_r = kwargs['triang_r']
        tpc = ax.tripcolor(triang_r, field, cmap=cmap, vmax = vmax, vmin = vmin, shading=shading)
    else:
        tpc = ax.tripcolor(triang, field, cmap=cmap, vmax = vmax, vmin = vmin,shading=shading)
        
    #cb =fig.colorbar(tpc)
    #cb.remove()
    #cb.ax.tick_params(labelsize=20)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    # Colorbar setting
    from mpl_toolkits.axes_grid1 import make_axes_locatable 
    divider = make_axes_locatable(ax)
    colax = divider.append_axes("right", size="1%", pad=0.2)#
    form = '%.1E'#'%1.2f'
    cbar = plt.colorbar(tpc, cax=colax, format=form, label=title) # 0.15 label=title,

    return tpc, ax

def plot_mode_fromfile(filename, varname, meshfile, **kwargs):
    """Load FELiCS fluctuation from file and plot it.

    Parameters
    ----------
    filename : string
        Full path of the solution file exported by FELiCS
    varname : string
        Name of the fluctuation variable to plot
    meshfile : string
        Full path of the mesh file exported by FELiCS
    **kwargs : dict
        The keyword arguments are passed to `plot_field()`

    Returns
    -------
    ax : matplotlib.axes
        ax object containing the plot
    tpc : matplotlib.collections.PolyQuadMesh
        pcolor object

    """
        
    if not 'quantity' in kwargs:
        quantity = 'magnitude'
    else:
        quantity = kwargs['quantity']

    # Get mesh connectivity and passing it to the 
    mesh_file_handle = h5py.File(meshfile, 'r')
    x = np.array(mesh_file_handle['coordinates/x'])
    y = np.array(mesh_file_handle['coordinates/y'])
    connectivity =  mesh_file_handle['cells']['triangles']
    triang = tri.Triangulation(x, y, triangles=connectivity)
    kwargs['triang'] = triang

    # Get mode data
    Out_file_handle1 = h5py.File(filename, 'r')
    if quantity in ['magnitude', 'angle']:
        varplot = np.array(Out_file_handle1['fluctuation']['0']['pointData'][varname[0]][quantity])
    elif quantity in ['real']:
        mag = np.array(Out_file_handle1['fluctuation']['0']['pointData'][varname[0]]['magnitude'])
        ang = np.array(Out_file_handle1['fluctuation']['0']['pointData'][varname[0]]['angle'])
        varplot = mag*np.cos(ang)
    elif quantity in ['imaginary']:
        mag = np.array(Out_file_handle1['fluctuation']['0']['pointData'][varname[0]]['magnitude'])
        ang = np.array(Out_file_handle1['fluctuation']['0']['pointData'][varname[0]]['angle'])
        varplot = mag*np.sin(ang)
    else:
        print("Attribute of mode shape not available.")
        exit()

    # Plot
    tpc, ax = plot_field(x, y, varplot, **kwargs)
    
    return
figure = plt.figure(figsize=(8, 6))
ax1 = plt.subplot(3,1,1)
plot_mode_fromfile(case_path+response_file,['ux'], case_path+mesh_file, quantity = 'real', ax=ax1, xlims=[-2,8], title='Real($u_x$)')
ax2 = plt.subplot(3,1,2)
plot_mode_fromfile(case_path+response_file,['ur'], case_path+mesh_file, quantity = 'real', ax=ax2, xlims=[-2,8], title='Real($u_r$)')
ax3 = plt.subplot(3,1,3)
plot_mode_fromfile(case_path+response_file,['p'], case_path+mesh_file, quantity = 'real', ax=ax3, xlims=[-2,8], title='Real($p$)')
plt.show()