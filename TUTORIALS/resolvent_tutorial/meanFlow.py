# Import required packages
import sys
import matplotlib.pyplot as plt
import numpy as np
import h5py

# Set the base path for data files

# File paths configuration
case = {
    'DATA_FILE': 'UMean.npy',  # Mean velocity data from LES data
    'MESH_FILE': 'FeliCS_mesh.npy',   # Mesh coordinates
    'SAVE_FILE': 'meanFlow.fel',           # Output file for FELiCS
}

# Physical parameters
physpar = {
    'Re': 8000,          # Reynolds number
    'L': 0.019/2,        # Reference length 
    'u_inf': 4*0.56842,  # Reference velocity 
}

# Load the mesh and base flow data
mesh = np.load(case['MESH_FILE'])
UMean = np.load(case['DATA_FILE'])

# Extract and ravel coordinates and velocities
x, y = mesh[:, :, 0].ravel(), mesh[:, :, 1].ravel()
u, v = UMean[:, :, 0].ravel(), UMean[:, :, 1].ravel()

# Normalize coordinates by reference length and velocities by reference velocity
x = x / physpar['L']
y = y / physpar['L']
ux = u / physpar['u_inf']
ur = v / physpar['u_inf']

nu_mol = (1 / physpar['Re']) * np.ones(len(x))
# NOTE: an eddy viscosity field can be defined as nuturb

# Create figure with the mean flow
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 6))
sc1 = ax1.scatter(x, y, c=ux)
ax1.set_title('Forcing Domain')
plt.colorbar(sc1, ax=ax1)
sc2 = ax2.scatter(x, y, c=ur)
ax2.set_title('Response Domain')
plt.colorbar(sc2, ax=ax2)

# Upstream sponge parameters
x0inlet = -2  # Start of inlet sponge region
a0 = 1 / (x[0] - x0inlet)**2  # Coefficient for quadratic sponge

# Create inlet sponge profile
spg_inlet = a0 * (x - x0inlet)**2
id_out = np.argwhere(x > x0inlet).ravel()  # Find points outside sponge region
spg_inlet[id_out] = 0  # Set sponge to zero outside defined region

# Downstream sponge parameters
x0outlet = 12  # Start of outlet sponge region
coeffOutlet = 1  # Coefficient scaling
a1 = coeffOutlet / (x[-1] - x0outlet)**2  # Coefficient for quadratic sponge

# Create outlet sponge profile
spg_outlet = a1 * (x - x0outlet)**2
id_out = np.argwhere(x < x0outlet).ravel()  # Find points outside sponge region
spg_outlet[id_out] = 0  # Set sponge to zero outside defined region

# Combine sponge regions
spg = spg_inlet + spg_outlet

# Plot the sponge regions
fig1, ax = plt.subplots()
tpc = plt.scatter(x, y, c=spg)
plt.axis('equal')
plt.colorbar(tpc)

# Define regions for forcing and response domains
Wforcing = np.where((x > -4) & (x < 12), 1, 0)
Wresponse = np.where((x > -4) & (x < 12), 1, 0)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 6))
sc1 = ax1.scatter(x, y, c=Wforcing)
ax1.set_title('Forcing Domain')
plt.colorbar(sc1, ax=ax1)
sc2 = ax2.scatter(x, y, c=Wresponse)
ax2.set_title('Response Domain')
plt.colorbar(sc2, ax=ax2)

# Create and save the FELiCS mean flow file
print(f'-- Saving data to file: {case["SAVE_FILE"]}')

hf = h5py.File(case['SAVE_FILE'], 'w')
# Save coordinates and velocity components
hf.create_dataset('/MeanFlow/x', data=x)
hf.create_dataset('/MeanFlow/r', data=y)
hf.create_dataset('/MeanFlow/ux', data=ux)
hf.create_dataset('/MeanFlow/ur', data=ur)

# Save viscosity field
hf.create_dataset('/MeanFlow/nulam', data=nu_mol)
#hf.create_dataset('/MeanFlow/nuturb', data=nuturb)
# Save sponge and domain information
hf.create_dataset('/MeanFlow/spg', data=spg)
hf.create_dataset('/MeanFlow/responseDomain', data=Wresponse)
hf.create_dataset('/MeanFlow/forcingDomain', data=Wforcing)

print('-- Data saved successfully!')