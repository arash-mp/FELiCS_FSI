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
import h5py
import numpy as np
import vtk
import meshio
import time
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
from scipy import interpolate

# Interpolation of large 3D meshes in FELiCS takes too long
# In this script the interpolation is done using linearNDinterpolator instead of griddata
# The mean fields are read from a vtk file but this can easily be changed to a hdf5 file
# I tried accelerating the interpolation by blocking the domain and interpolate in parallel for each block. 
# There are still errors in the results when using this method which is why it is commented out.
# Instead linearND allows to do the triangulation only once and then interpolating each field in parallel (done by using ThreadPoolExecutor)
# Masking is commented out but can be used to reduce the size of the domain similar to a bounding box
# Data size can further be reduced by skipping points in the domain (n_s)
# the "normalised" boolean is used to determine wether the domain is to be normalised or not

########################################################################################################
########################################## Read vtk file ###############################################
########################################################################################################
start_time_reading = time.time()

# set up reader
filename = 'CC_data.vtk'
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(filename) # file name of the vtk file
reader.Update()

print('Reading data from vtk file: ',filename,'...')

data = reader.GetOutput() # get the output

points = data.GetPoints() # retrieve points of the mesh
point_coords = np.zeros((points.GetNumberOfPoints(), 3)) # for 3D

for i in range(points.GetNumberOfPoints()):
    point_coords[i, :] = points.GetPoint(i)

print("Read the fields...")
point_data = data.GetPointData()
velocity_array = point_data.GetArray("Velocity")

npoints = velocity_array.GetNumberOfTuples() # number of points
ncomp = velocity_array.GetNumberOfComponents() # number of components

velocity_data = np.zeros((npoints, ncomp)) 
for i in range(npoints):
    velocity_data[i, :] = velocity_array.GetTuple3(i) # fill velocity array with data from the vtk

vis_array = point_data.GetArray("Vis")
vis_data = np.zeros((npoints, 1))
for i in range(npoints):
    vis_data[i, :] = vis_array.GetTuple1(i) # fill viscosity array with data from the vtk

q_dot_array = point_data.GetArray("HeatReleaseRate")
q_dot_data = np.zeros((npoints, 1))
for i in range(npoints):
    q_dot_data[i, :] = q_dot_array.GetTuple1(i) # fill heat release rate array with data from the vtk   

# phi_forcing_r = np.zeros((npoints, 1))
# phi_forcing_r[:, 0] = 1
phi_forcing_r = np.ones((npoints, 1))

density_array = point_data.GetArray("density") # fill density array with data from the vtk
density_data = np.zeros((npoints, 1))
for i in range(npoints):
    density_data[i, :] = density_array.GetTuple1(i)

# Extract boundary nodes
surface_filter = vtk.vtkDataSetSurfaceFilter()
surface_filter.SetInputData(data)
surface_filter.Update()

boundary_mesh = surface_filter.GetOutput()
boundary_points = set()
boundary_points_np = np.array([boundary_mesh.GetPoint(i) for i in range(boundary_mesh.GetNumberOfPoints())])
boundary_points = set(map(tuple, boundary_points_np))
boundary_indices = [i for i in range(points.GetNumberOfPoints()) if tuple(points.GetPoint(i)) in boundary_points]

all_points = set(tuple(points.GetPoint(i)) for i in range(points.GetNumberOfPoints()))
internal_points = all_points - boundary_points
internal_indices = [i for i in range(points.GetNumberOfPoints()) if tuple(points.GetPoint(i)) in internal_points]

end_time_reading = time.time()
time_reading = end_time_reading - start_time_reading
print(f'Reading from the vtk file and writing to arrays took: {time_reading:.4f} seconds.')
########################################################################################################
######################################### Read .msh file ###############################################
########################################################################################################
print('Reading the mesh file...')
normalised = True
mesh = meshio.read(f"final_tf.msh")
if normalised:
    x_mm = mesh.points[:, 0]/np.max(mesh.points[:, 0])  # x-coordinates
    y_mm = mesh.points[:, 1]/np.max(mesh.points[:, 0])  # y-coordinates
    z_mm = mesh.points[:, 2]/np.max(mesh.points[:, 0])  # z-coordinates
    print(f"Mesh loaded with {x_mm.shape[0]} nodes.")
else:
    x_mm = mesh.points[:, 0]  # x-coordinates
    y_mm = mesh.points[:, 1]  # y-coordinates
    z_mm = mesh.points[:, 2]  # z-coordinates

print(f"Mesh loaded with {x_mm.shape[0]} nodes.")

########################################################################################################
####################################### Naming of variables ############################################
########################################################################################################
print('Naming variables...')
if normalised:
    x_f = point_coords[:, 0].reshape(npoints,1)/np.max(mesh.points[:, 0])
    y_f = point_coords[:, 1].reshape(npoints,1)/np.max(mesh.points[:, 0])
    z_f = point_coords[:, 2].reshape(npoints,1)/np.max(mesh.points[:, 0])
else:
    x_f = point_coords[:, 0].reshape(npoints,1)
    y_f = point_coords[:, 1].reshape(npoints,1)
    z_f = point_coords[:, 2].reshape(npoints,1)
ux_f = velocity_data[:,0].reshape(npoints,1)
uy_f = velocity_data[:,1].reshape(npoints,1)
uz_f = velocity_data[:,2].reshape(npoints,1)
nuturb_f = vis_data + np.max(vis_data)/80
D_phi_f = nuturb_f/0.9
phi_forcing_r_f = phi_forcing_r
rho_f = density_data
q_dot_f = q_dot_data

########################################################################################################
###################################### Reducing data density ###########################################
########################################################################################################

n_s = 5
print(f'Reducing data density by factor {n_s}...')
x = np.concatenate((x_f[boundary_indices], x_f[internal_indices][::n_s]))
y = np.concatenate((y_f[boundary_indices], y_f[internal_indices][::n_s]))
z = np.concatenate((z_f[boundary_indices], z_f[internal_indices][::n_s]))
ux = np.concatenate((ux_f[boundary_indices], ux_f[internal_indices][::n_s]))
uy = np.concatenate((uy_f[boundary_indices], uy_f[internal_indices][::n_s]))
uz = np.concatenate((uz_f[boundary_indices], uz_f[internal_indices][::n_s]))
D_phi = np.concatenate((D_phi_f[boundary_indices], D_phi_f[internal_indices][::n_s]))
phi_forcing_r = np.concatenate((phi_forcing_r_f[boundary_indices], phi_forcing_r_f[internal_indices][::n_s]))
nuturb = np.concatenate((nuturb_f[boundary_indices], nuturb_f[internal_indices][::n_s]))
rho = np.concatenate((rho_f[boundary_indices], rho_f[internal_indices][::n_s]))
q_dot = np.concatenate((q_dot_f[boundary_indices], q_dot_f[internal_indices][::n_s]))

x_m = x_mm
y_m = y_mm
z_m = z_mm
########################################################################################################
######################################## Masking meanfields ############################################
########################################################################################################
# print('Masking of the mean field...')
# z_max = 0.264
# z_min = 0.169

# # Filter based on z condition
# mask = ((z_f >= z_min) & (z_f <= z_max))
# mask_m = ((z_mm >= z_min) & (z_mm <= z_max))
# # Masking 
# x = x_f[mask]
# y = y_f[mask]
# z = z_f[mask]
# ux = ux_f[mask]
# uy = uy_f[mask]
# uz = uz_f[mask]
# D_phi = D_phi_f[mask]
# phi_forcing_r = phi_forcing_r_f[mask]
# nuturb = nuturb_f[mask]
# rho = rho_f[mask]
# q_dot = q_dot_f[mask] 
# x_m = x_mm[mask_m]
# y_m = y_mm[mask_m]
# z_m = z_mm[mask_m]


########################################################################################################
###################################### Blocking of meanfield ###########################################
########################################################################################################
# print('Blocking the mean fields...')
# start_time_blocking = time.time()

# x_blocking = 3 # divisions in x direction
# z_b = np.min(z) + (np.max(z) - np.min(z)) / 2 # find midpoint in z
# x_bn = (np.max(x) - np.min(x)) / x_blocking # defining width of blocks in x direction

# # Create an array of x-block boundaries
# x_b = np.linspace(np.min(x), np.min(x) + (x_blocking - 1) * x_bn, x_blocking)

# overlap = 0.05 # percentage of overlap

# i_x = [np.where((x > x_b[i] - (overlap*x_bn)) & (x <= x_b[i + 1] + (overlap*x_bn)) if i < x_blocking - 1 else x > x_b[i] - (overlap*x_bn))[0]
#        for i in range(x_blocking)]
# i_x_m = [np.where((x_m > x_b[i] ) & (x_m <= x_b[i + 1]) if i < x_blocking - 1 else x_m > x_b[i])[0]
#        for i in range(x_blocking)]

# # identify indices for y and z conditions
# i_yn = np.where(y <= 0 + (overlap*np.max(y)))[0]
# i_yp = np.where(y >= 0 - (overlap*np.max(y)))[0]
# i_y = [i_yn, i_yp]


# i_zn = np.where(z <= z_b + (overlap*((np.max(z) - np.min(z)) / 2)))[0]
# i_zp = np.where(z >= z_b - (overlap*((np.max(z) - np.min(z)) / 2)))[0]
# i_z = [i_zn, i_zp]

# # for the target mesh no overlap
# i_yn_m = np.where(y_m <= 0)[0] 
# i_yp_m = np.where(y_m > 0)[0]
# i_y_m = [i_yn_m, i_yp_m]

# i_zn_m = np.where(z_m <= z_b)[0]
# i_zp_m = np.where(z_m > z_b)[0]
# i_z_m = [i_zn_m, i_zp_m]

# # array of strings for efficient handling of fields in loop
# fields = ["x", "y", "z", "ux", "uy", "uz", "D_phi", "phi_forcing_r", "nuturb", "rho", "q_dot"]
# fields_m = ["x_m", "y_m", "z_m"]

# # initialize dictionary for the blocked data
# data_blocks = {field: [] for field in fields}
# data_blocks_m = {field: [] for field in fields_m}

# # loop through each block
# n_blocks = 0
# for i in range(len(i_x)):
#     for j in range(len(i_y)):
#         for k in range(len(i_z)):
#             # indice for current block
#             i_c = np.intersect1d(np.intersect1d(i_x[i], i_y[j]), i_z[k])
#             i_c_m = np.intersect1d(np.intersect1d(i_x_m[i], i_y_m[j]), i_z_m[k])
#             if i_c.size > 1:
#                 # store data for each block
#                 for field in fields:
#                     data_blocks[field].append(eval(field)[i_c]) # evaluates each field at the indices at the current block
#                 for field_m in fields_m:
#                     data_blocks_m[field_m].append(eval(field_m)[i_c_m]) # for target mesh
#                 n_blocks += 1
#                 print(f"Block ({n_blocks}) covers {len(i_c_m)} data points.")

# # Check if total number of points across all blocks equals the original size
# total_points = sum([len(np.intersect1d(np.intersect1d(i_x[i], i_y[j]), i_z[k])) for i in range(len(i_x)) for j in range(len(i_y)) for k in range(len(i_z))])
# print(f"Total points covered by all blocks: {total_points}")
# total_points_m = sum([len(np.intersect1d(np.intersect1d(i_x_m[i], i_y_m[j]), i_z_m[k])) for i in range(len(i_x_m)) for j in range(len(i_y_m)) for k in range(len(i_z_m))])
# print(f"Total points covered by all blocks: {total_points_m}")

# end_time_blocking = time.time()
# time_blocking = end_time_blocking - start_time_blocking
# print(f"Data  successfully blocked into {n_blocks} blocks in {time_blocking:.4f} seconds.")
# input('wait')
########################################################################################################
#################################### Parallelise interpolation #########################################
########################################################################################################
# print('Interpolating mean fields to mesh...')

# def interpolate_block(i):
#     start_time_tria = time.time()
    
#     # Gather data for the current block
#     x_block = data_blocks['x'][i]
#     y_block = data_blocks['y'][i]
#     z_block = data_blocks['z'][i]
#     points_block = np.column_stack((x_block, y_block, z_block))
    
#     # Mean fields to interpolate
#     mean_fields = np.column_stack((
#         data_blocks['ux'][i], data_blocks['uy'][i], data_blocks['uz'][i],
#         data_blocks['rho'][i], data_blocks['phi_forcing_r'][i],
#         data_blocks['nuturb'][i], data_blocks['q_dot'][i], data_blocks['D_phi'][i]
#     ))
    
#     # Target points
#     target_points = np.column_stack((
#         data_blocks_m['x_m'][i], data_blocks_m['y_m'][i], data_blocks_m['z_m'][i]
#     ))
    
#     # Perform triangulation and interpolation
#     triangulation = Delaunay(points_block)
#     interpolator = LinearNDInterpolator(triangulation, mean_fields)
#     interpolated_values = interpolator(target_points)
    
#     end_time_tria = time.time()
#     print(f'Interpolating block {i + 1} took: {end_time_tria - start_time_tria:.4f} seconds.')
    
#     return interpolated_values

# # Start the parallel processing
# start_time_interpolation = time.time()
# interpolated_fields_all_blocks = []

# with ProcessPoolExecutor() as executor:
#     results = list(executor.map(interpolate_block, range(n_blocks)))

# # Collect results
# interpolated_fields_all_blocks = list(results)
# for i, block in enumerate(interpolated_fields_all_blocks):
#     print(f"Block {i+1} shape: {block.shape}")
# end_time_interpolation = time.time()
# time_interpolation = end_time_interpolation - start_time_interpolation
# print(f'Overall interpolation of all blocks took: {time_interpolation:.4f} seconds with {x.shape} nodes interpolated on {x_m.shape} nodes of the target mesh.')

# total_points = sum([block.shape[0] for block in interpolated_fields_all_blocks])
# print(f"Total interpolated points: {total_points}")

########################################################################################################
################################# Interpolating to mesh (blocked) ######################################
########################################################################################################

print('Triangulating...')
points_xyz = np.column_stack((x,y,z)) # array of coordinates for the vtk mesh
points_xyz_m = np.column_stack((x_m,y_m,z_m)) # array of coordinates for the target mesh

start_tria = time.time()
triangulation = Delaunay(points_xyz) # triangulate the vtk mesh once
end_tria = time.time()
time_tria = end_tria - start_tria
print(f'Triangulating took: {time_tria:.4f} seconds.')

def interpolate_with_fallback(points_xyz, field_data, mesh_points, triangulation, field_name):
    # function to interpolate the fields with fallback for NaN values
    print(f"Interpolating {field_name} with shape {field_data.shape}...")
    interpolator = LinearNDInterpolator(triangulation, field_data)
    interpolated = interpolator(mesh_points)
    
    if np.any(np.isnan(interpolated)):
        print(f"Handling NaN values for {field_name} with griddata...")
        fallback = interpolate.griddata(points_xyz, field_data, mesh_points, method='nearest') # NearestNDInterpolator maybe faster (not tested yet)
        interpolated[np.isnan(interpolated)] = fallback[np.isnan(interpolated)] # replace NaN values with nearest neighbour
    return interpolated

interpolated_fields = {} # list of all the fields to interpolate in parallel
fields_to_interpolate = {
    "ux": ux,
    "uy": uy,
    "uz": uz,
    "nuturb": nuturb,
    "q_dot": q_dot,
    "rho": rho
} 

for field_name, field_data in fields_to_interpolate.items():
    if len(points_xyz) != len(field_data):
        raise ValueError(
            f"Mismatch in lengths: point_coords has {len(points_xyz)} points "
            f"but {field_name} has {len(field_data)} values."
        )

print("Interpolating fields onto target mesh using VTK triangulation...")
start_inter = time.time()

interpolated_fields = {}
interpolated_fields_griddata = {}

with ThreadPoolExecutor(max_workers=len(fields_to_interpolate)) as executor:
    results = executor.map(
        lambda args: interpolate_with_fallback(*args),
        zip([points_xyz] * len(fields_to_interpolate),
            fields_to_interpolate.values(),
            [points_xyz_m] * len(fields_to_interpolate),
            [triangulation] * len(fields_to_interpolate),
            fields_to_interpolate.keys())
    )
    interpolated_fields = {field_name: result for field_name, result in zip(fields_to_interpolate.keys(), results)}

end_inter = time.time()
time_inter = end_inter - start_inter
print(f'Interpolating all fields took: {time_inter:.4f} seconds.')

########################################################################################################
####################################### Write the .fel file ############################################
########################################################################################################
# input('wait')
print('Writing data to .fel file...')
# FELiCS file name
felFile = "meanflow_final_tf_msh.fel"

# Save mean field data to a HDF5  .fel file
with h5py.File(felFile, 'w') as f:
    f.create_dataset('/MeanFlow/x', data=x_m)
    f.create_dataset('/MeanFlow/y', data=y_m)
    f.create_dataset('/MeanFlow/z', data=z_m)
    f.create_dataset('/MeanFlow/D_phi', data=interpolated_fields.get("nuturb")/0.9)
    f.create_dataset('/MeanFlow/phi_forcing_r', data=np.ones_like(x_m))
    # for field, data in field_data.items():
    for field, data in interpolated_fields.items():    
        dataset_path = f'/MeanFlow/{field}'
        f.create_dataset(dataset_path, data=data)

print('Done.')
