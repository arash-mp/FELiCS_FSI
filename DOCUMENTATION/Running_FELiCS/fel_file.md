# fel file in FELiCS
In progress 

Explains how to write the .fel file containing the mase fglow information.

''' python
import  h5py
import  numpy as np

case = {
    'PATH_MEANFLOW':    'myFELiCSdir/', 
    'DATA_FILE':        'input_base_flow/', #input base flow data
    'MESH_FILE':        'input_mesh', # input mesh
    'SAVE_FILE':        'myFELiCScase/Base_flow/base_flow.fel', #saving directory
}

hf = h5py.File(save_dir, 'w')
hf.create_dataset('/MeanFlow/x', data=x)
hf.create_dataset('/MeanFlow/r', data=y)
hf.create_dataset('/MeanFlow/ux', data=ux)
hf.create_dataset('/MeanFlow/ur', data=ur)
hf.create_dataset('/MeanFlow/ut', data=ut)
hf.create_dataset('/MeanFlow/nulam', data=nu_mol)

# One can include here an eddy viscosity field
hf.create_dataset('/MeanFlow/nuturb', data=nut)

# One can include here a sponge field
hf.create_dataset('/MeanFlow/spg', data=spg)
    
hf.create_dataset('/MeanFlow/responseDomain', data=Wresponse.restrictor)
    hf.create_dataset('/MeanFlow/forcingDomain', data=Wforcing.restrictor)
hf.close()
'''