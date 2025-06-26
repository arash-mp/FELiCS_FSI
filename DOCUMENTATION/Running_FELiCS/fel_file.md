# fel file in FELiCS
FELiCS prcesses the mean flow information with a .fel file using the .h5 format. It contains the mesh, mean flow field as well as the forcing and response domains.
It is also possible to add a sponge field and an eddy viscosity field. 

Here is a python script example of how the fel file can be written: 
``` python
import  h5py

case = {
    'PATH_MEANFLOW':    'myFELiCSdir/', 
    'DATA_FILE':        'input_base_flow/', #input base flow data
    'MESH_FILE':        'input_mesh', # input mesh
    'SAVE_FILE':        'myFELiCScase/Base_flow/base_flow.fel', #saving directory
}

# The data from the mesh and base flow are loaded

hf = h5py.File(case['SAVE_FILE'], 'w')
hf.create_dataset('/MeanFlow/x', data=x)
hf.create_dataset('/MeanFlow/r', data=y)
hf.create_dataset('/MeanFlow/ux', data=ux)
hf.create_dataset('/MeanFlow/uy', data=uy)
hf.create_dataset('/MeanFlow/uz', data=uz)
hf.create_dataset('/MeanFlow/nulam', data=nu_mol)

# One can include here an eddy viscosity field
hf.create_dataset('/MeanFlow/nuturb', data=nut)

# One can include here a sponge field
hf.create_dataset('/MeanFlow/spg', data=spg)
    
hf.create_dataset('/MeanFlow/responseDomain', data=Wresponse)
    hf.create_dataset('/MeanFlow/forcingDomain', data=Wforcing)
hf.close()
```