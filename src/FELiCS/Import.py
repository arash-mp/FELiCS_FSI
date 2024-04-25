#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
# **********************************************************************
'''
# **********************************************************************
# * This file deals with reading flow data from various sources and
# * Interpolating it on the computational grid.
# * This file was created by Thomas L. Kaiser. Significant contributions
# * were done by...
# **********************************************************************
'''

import numpy as np
def ExpandForAverage(coordinates,param):
    ''' This function extends the fenics grid to three dimensions to prepare for interpolation in 3D grids. Afet interpolation use the function ContractAfterAverage to project the 3D data again on the 2D FELiCS mesh.
        \t Input:
        \t\t -coordinates: Coordinates of the FELiCS mesh
        \t\t -param: FELiCS parameter object
        \t Ouput:
        \t\t -coordinates_out: FELiCS mesh coordinates expanded to 3 dimensions
        '''
    n_cuts=10
    size=np.shape(coordinates)[0]
    b=np.zeros((size,3))
    b[:,0]=coordinates[:,0]
    b[:,1]=coordinates[:,1]
    coordinates=b

    temp_coordinates=np.zeros(np.shape(coordinates))
    coordinates_out=np.zeros((0,3))
    if param.Case.CoordinateSystem=='Cylindrical':
        for angle in np.linspace(2*np.pi/n_cuts,2*np.pi,n_cuts):
            temp_coordinates[:,0]=coordinates[:,0]
            temp_coordinates[:,1]=coordinates[:,1]*np.cos(angle)
            temp_coordinates[:,2]=coordinates[:,1]*np.sin(angle)
            coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
    elif param.Case.CoordinateSystem=='Cartesian':
        zmin=-0.0025
        zmax=0.0025
        for z in np.linspace(zmin,zmax,n_cuts):
            temp_coordinates[:,0]=coordinates[:,0]
            temp_coordinates[:,1]=coordinates[:,1]
            temp_coordinates[:,2]=coordinates[:,1]*0+z
            coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
    return coordinates_out

def ContractAfterAverage(coordinates,param,vec_in):
    ## This function contracts the expanded interpolated data from the expaned grid (see function ExpandForZaimuthalAverage) back
    ## to the 2D base grid
    n_cuts=10
    length=np.shape(coordinates)[0]
    vec_out=np.zeros(length)
    for i in range(n_cuts):
        vec_out=vec_out+vec_in[i*length:i*length+length]
    vec_out=vec_out/n_cuts
    return vec_out


def importHDF5File(param,FEMSpaces):
    from fenics import HDF5File,Function
    ## Define Dictionaries
    MeanFlowDict={}
    notInFileList=[]
    ## Get mesh data
    mesh=FEMSpaces.P2.mesh()
    # If hdf5 file is used (so far implemented for fenics constructed hdf5 files)
    hdf5file = HDF5File(mesh.mpi_comm(),param.FlowInput.MeanFlowFilePath, 'r')
    for name in param.Case.getMeanFlowFieldNames():
        if name[0] == 'u':
            MeanFlowDict[name]=Function(FEMSpaces.FunctionSpaceVectorVelocity)
        else:
            MeanFlowDict[name]=Function(FEMSpaces.P2)
        hdf5file.read(MeanFlowDict[name], name)
    print(MeanFlowDict.keys())
    return MeanFlowDict
