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
# Standard libraries
from functools import partial
import multiprocessing
import os
import pickle
import subprocess
import sys
import time

# Third party libraries
from dolfinx.fem import Function
import numpy as np
from scipy.interpolate import griddata

def get_last_git_commit():
    """
    Get the latest git commit hash of the repository.

    Returns
    -------
    str
        The most recent commit hash as a string.
    """
    FELiCSPathname = os.path.dirname(sys.argv[0])

    commit = subprocess.Popen(
    ['git','--git-dir',FELiCSPathname+'/../.git', 'rev-parse', 'HEAD'],
    shell=False,
    stdout=subprocess.PIPE,
    )
    commit = commit.communicate()[0].strip().decode('ascii')
    return commit


#####################################
#### Deprecated code ################
#####################################
# To revive, create a working test for it.


#def getSpeciesListSolution(param):
#    """
#    Get the list of species to be solved for, based on the chemistry model.
#
#    Parameters
#    ----------
#    param : object
#        Object containing the chemistry model and optional additional species.
#
#    Returns
#    -------
#    list of str
#        List of species names included in the solution.
#    """
#    #import parameters as param
#    SpeciesList=[]
#    if param.ChemistryModel=='OneStep':
#        SpeciesList.append('O2' )
#        SpeciesList.append('CH4' )
#    elif param.ChemistryModel=='BFER':
#        SpeciesList.append('CH4' )
#        SpeciesList.append('O2' )
#        SpeciesList.append('CO2' )
#        SpeciesList.append('CO')
#        SpeciesList.append('H2O')
#    # Add user specified additional species
#    for specie in param.additionalSpecies:
#        SpeciesList.append(specie)
#    return SpeciesList
#
#def getSpeciesListMean(param):
#    """
#    Get the list of species required in the mean fields based on the chemistry model.
#
#    Parameters
#    ----------
#    param : object
#        Object containing the chemistry model and optional additional species.
#
#    Returns
#    -------
#    list of str
#        List of species names included in the mean fields.
#    """
#    #import parameters as param
#    SpeciesList=[]
#    if param.ChemistryModel=='OneStep':
#        #SpeciesList.append(param.ReactionProgressMarker )
#        SpeciesList.append('CH4' )
#        SpeciesList.append('O2')
#        SpeciesList.append('CO2' )
#        SpeciesList.append('CO')
#    elif param.ChemistryModel=='BFER':
#        SpeciesList.append('CH4' )
#        SpeciesList.append('O2' )
#        SpeciesList.append('CO2' )
#        SpeciesList.append('CO')
#        SpeciesList.append('H2O')
#        SpeciesList.append('N2')
#    # Add user specified additional species
#    for specie in param.additionalSpecies:
#        SpeciesList.append(specie)
#    return SpeciesList
#
#def getReactionList(param):
#    """
#    Get the list of chemical reactions defined for the selected chemistry model.
#
#    Parameters
#    ----------
#    param : object
#        Contains the name of the chemistry model.
#
#    Returns
#    -------
#    list of dict
#        Each dictionary contains details about one chemical reaction.
#    """
#    #import parameters as param
#    List=[]
#    if param.ChemistryModel=='OneStep':
#        tempDict={}
#        tempDict['educts']=['O2', 'CH4']
#        tempDict['educt_mole_factor']=[2,1]
#        tempDict['products']=['H2O','CO2']
#        tempDict['product_mole_factor']=[2,1]
#        tempDict['A']=6.7e9
#        tempDict['E_a']=48400
#        tempDict['R']=1.987
#        tempDict['a']=1.3
#        tempDict['b']=0.2
#        tempDict['beta']=0.0
#        tempDict['h0']=-800e3
#        tempDict['modelname']='OneStep'
#        List.append(tempDict)
#    elif param.ChemistryModel=='BFER':
#        #See notations and parameters in B.Franzelli et al. Combustion and flame. (2012)
#
#        #Reation O2CH4 : CH4+1.5O2 => CO + 2H2O
#        tempDict={}
#        tempDict['educts']=['O2', 'CH4']
#        tempDict['educt_mole_factor']=[1.5,1.0]
#        tempDict['products']=['H2O','CO']
#        tempDict['product_mole_factor']=[2.0,1.0]
##       tempDict['A']=4.9e9  #cgs
#        tempDict['A']=4.9e9*(1/10**0.9)  #unit conversion from cgs to SI.
#        tempDict['beta']=0.0
#        tempDict['E_a']=3.55e4
#        tempDict['R']=1.987
#        #correction function f1
#        tempDict['phi0']=1.1
#        tempDict['segma0']=0.09
#        tempDict['B']=0.37
#        tempDict['phi1']=1.13
#        tempDict['segma1']=0.03
#        tempDict['C']=6.7
#        tempDict['phi2']=1.6
#        tempDict['segma2']=0.22
#        tempDict['b']=0.5 #nC4
#        tempDict['a']=0.65 #nO2_1
#        tempDict['h0']=-519.26e3
#        tempDict['modelname']='BFER'
#        List.append(tempDict)
#
#
#        #Reaction O2CO : CO+0.5O2 => CO2
#        tempDict={}
#        tempDict['educts']=['O2', 'CO']
#        tempDict['educt_mole_factor']=[0.5,1.0]
#        tempDict['products']=['CO2']
#        tempDict['product_mole_factor']=[1.0]
##       tempDict['A']=2.0e8 #cgs
#        tempDict['A']=2.0e5  #unit conversion from cgs to SI.
#        tempDict['beta']=0.7
#        tempDict['E_a']=1.2e4
#        tempDict['R']=1.987
#        #correction function f2
#        tempDict['phi0']=0.95
#        tempDict['segma0']=0.08
#        tempDict['B']=2.5e-5
#        tempDict['phi1']=1.3
#        tempDict['segma1']=0.04
#        tempDict['C']=0.0087
#        tempDict['phi2']=1.2
#        tempDict['segma2']=0.04
#        tempDict['phi3']=1.2
#        tempDict['segma3']=0.05
#        tempDict['b']=1.0 #nCO
#        tempDict['a']=0.5 #nO2_2
#        tempDict['h0']=-282.98e3
#        tempDict['modelname']='BFER'
#        List.append(tempDict)
#    return List
#
#def getReverseReactionList(param):
#    """
#    Get a list of reverse reactions for the specified chemistry model.
#
#    Parameters
#    ----------
#    param : object
#        Object containing chemistry model information.
#
#    Returns
#    -------
#    list of dict
#        List of reverse chemical reactions.
#    """
#    List=[]
#    if param.ChemistryModel=='OneStep':
#        pass
#    elif param.ChemistryModel=='BFER':
#        #Reaction O2CO : CO2 => 0.5O2+CO
#        tempDict={}
#        tempDict['educts']=['CO2']
#        tempDict['educt_mole_factor']=[1.0]
#        tempDict['products']=['O2','CO']
#        tempDict['product_mole_factor']=[0.5,1]
#        tempDict['fa']=0.5 #nO2_2
#        tempDict['fb']=1.0 #nCO
#        tempDict['rc']=1.0 #nCO2
#        tempDict['h0']=519.26e3
#        tempDict['modelname']='BFER'
#
#        List.append(tempDict)
#
#    return List
#
#def getMolecularMass():
#    """
#    Get a dictionary of molecular masses for key chemical species.
#
#    Returns
#    -------
#    dict
#        Mapping of species names to their molecular masses.
#    """
#    MolecularMass={}
#    MolecularMass['O2']=0.0319988e0
#    MolecularMass['CH4']=0.0160423e0
#    MolecularMass['H2O']=3.002628e-02
#    MolecularMass['CO']=0.0280140e0
#    MolecularMass['N2']=0.0280134e0
#    MolecularMass['CO2']=0.0440098e0
#
#    return MolecularMass
#
#def getReactionName(Reaction):
#    """
#    Get the name of a reaction by concatenating its educt species.
#
#    Parameters
#    ----------
#    Reaction : dict
#        Dictionary containing a key 'educts' with a list of species names.
#
#    Returns
#    -------
#    str
#        Concatenated name of the reaction.
#    """
#    ReactionName=''
#    for specie in Reaction['educts']:
#        ReactionName=ReactionName+specie
#    return ReactionName
#
#def getJANAFTable():
#    """
#    Get thermodynamic data from the JANAF table for key species.
#
#    Returns
#    -------
#    dict
#        Dictionary of JANAF polynomial coefficients for multiple species.
#    """
#    #Ref temperature T=0 K
#    #cp=R*(a1+a2*T+a3*T**2+a4*T**3+a5*T**4)
#    #h=R*(a1*T+a2*T**2/2+a3*T**3/3+a4*T**4/4+a5*T**5/5+a6)
#    #s=R*(a1*ln(T)+a2*T+a3*T**2/2+a4*T**3/3+a5*T**4/4+a7)
#    #for T>1000K [a1, a2, a3, a4, a5, a6\
#    #for T<1000K a1, a2, a3, a4, a5, a6]
#    tempDict={}
#    tempDict['O2']=[3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14,-1.08845772E+03,5.45323129E+00,\
#         3.78245636E+00,-2.99673416E-03, 9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00]
#    tempDict['H2O']=[3.03399249E+00, 2.17691804E-03,-1.64072518E-07,-9.70419870E-11, 1.68200992E-14,-3.00042971E+04,4.96677010E+00,\
#         4.19864056E+00,-2.03643410E-03, 6.52040211E-06,-5.48797062E-09, 1.77197817E-12,-3.02937267E+04,-8.49032208E-01]
#    tempDict['CH4']=[7.48514950E-02, 1.33909467E-02,-5.73285809E-06, 1.22292535E-09,-1.01815230E-13,-9.46834459E+03,\
#         1.84373180E+01,5.14987613E+00,-1.36709788E-02, 4.91800599E-05,-4.84743026E-08, 1.66693956E-11,-1.02466476E+04,-4.64130376E+00 ]
#    tempDict['CO']=[2.71518561E+00, 2.06252743E-03,-9.98825771E-07, 2.30053008E-10,-2.03647716E-14,-1.41518724E+04,7.81868772E+00,\
#         3.57953347E+00,-6.10353680E-04, 1.01681433E-06,9.07005884E-10,-9.04424499E-13,-1.43440860E+04,3.50840928E+00]
#    tempDict['CO2']=[3.85746029E+00, 4.41437026E-03,-2.21481404E-06, 5.23490188E-10,-4.72084164E-14,-4.87591660E+04, 2.27163806E+00,\
#         2.35677352E+00, 8.98459677E-03,-7.12356269E-06,2.45919022E-09,-1.43699548E-13,-4.83719697E+04,9.90105222E+00 ]
#    return tempDict
#
#def index_2d(myList, v):
#    """
#    Find the 2D index of a value in a nested list.
#
#    Parameters
#    ----------
#    myList : list of list
#        The 2D list to search in.
#    v : object
#        The value to search for.
#
#    Returns
#    -------
#    tuple of int
#        Indices of the value in the list.
#    """
#    for i, x in enumerate(myList):
#        if v in x:
#            return (i, x.index(v))
#
#def rotateTheta(coordinates,Theta):
#    """
#    Rotate 3D coordinates in the YZ-plane by a given angle Theta.
#
#    This function takes coordinates (3D required) and calculates new coordinates by rotating the original
#    coordinates with respect to the angle theta. The first component is not altered; the second and third
#    are computed.
#
#    Parameters
#    ----------
#    coordinates : np.ndarray
#        Array of shape (n, 3) representing 3D points.
#    Theta : float
#        Angle by which to rotate the YZ components.
#
#    Returns
#    -------
#    np.ndarray
#        Rotated 3D coordinates.
#    """
#    # bring coordinates in right shape
#    size=np.shape(coordinates)[0]
#    b=np.zeros((size,3))
#    b[:,0]=coordinates[:,0]
#    b[:,1]=coordinates[:,1]
#    coordinates=b
#
#    temp_coordinates=np.zeros(np.shape(coordinates))
#    coordinates_out=np.zeros((0,3))
#    temp_coordinates[:,0]=coordinates[:,0]
#    temp_coordinates[:,1]=coordinates[:,1]*np.cos(Theta)
#    temp_coordinates[:,2]=coordinates[:,1]*np.sin(Theta)
#    coordinates_out=np.concatenate((coordinates_out,temp_coordinates))
#    return coordinates_out
#
#def ExtrudeFelicsGridToVTK(coords2D, Tri2D, angularSteps):
#    """
#    Convert a 2D grid into a 3D cylindrical grid by azimuthal extrusion.
#
#    This function expands a planar, 2D grid into the azimuthal direction. It requires the coordinates of
#    the cells' vertices, the number of planes in azimuthal direction (angularSteps), and the 2D grid triangulation.
#
#    Parameters
#    ----------
#    coords2D : np.ndarray
#        Coordinates of the 2D grid vertices.
#    Tri2D : np.ndarray
#        2D triangulation information (triangles).
#    angularSteps : int
#        Number of angular divisions in the azimuthal direction.
#
#    Returns
#    -------
#    tuple
#        Triangulation and coordinates for the resulting 3D grid.
#    """
#    ''' 
#    '''
#    # calculate coordinates of first new plane
#    coordsRotated=rotateTheta(coords2D,2*3.141/angularSteps)
#
#    #   bring coordinates into suiteable list format xyz2(Vertex_j) = [x(j), y(j), z(j)]
#    xyz2=[[k[0],k[1],k[0]*0] for k in coords2D]
#    xyz3D=[[k[0],k[1],k[2]] for k in coordsRotated]
#
#    # change triangle triangulation (3 indices per element) into triangulation for the first slice
#    # one slice is the triangulation between to flat, neighbouring planes
#    # triangulation is built by wedges (6 indices per element) <- Tri3
#    # triangluation indices of wedge element only needs to be extended by the same indices from triangles....
#    # with an addtion of the vertices in the in the 2D plane
#    Tri3 = []
#    Ntriangles = Tri2D.shape[0] # counting the amount of elements in 2D (amount if triangles)
#    Npoints = coords2D.shape[0] # counting the amount of vertices in 2D
#    for k in range(0,Ntriangles): # loop over all elements/triangles
#        Tri3.append(list(Tri2D[k])) # Tri[k] = [1,2,3]
#        Tri3[k].extend([kk+Npoints for kk in Tri2D[k]]) # Tri[k] = [1,2,3, 1+Npoints, 2+Npoints, 3+Npoints]
#
#    # coordinates corresponding to the vertices of first slice
#    # extension of 2D-ccordinates by earlier rotated coordinates
#    coords3D = xyz2+xyz3D
#
#    # Init. of Theta (counting up by dTheta every loop cycle) AND dTheta
#    Theta = dTheta = 2*np.pi/(angularSteps)
#    # loop (NOT starting from 0, since first "loop" was finished above) over
#    # remaining planes
#    for m in range(1,angularSteps):
#        Theta += dTheta
#
#        # similar as above: 2D-coordinates of origin-plane, only rotated by another Theta here!
#        coordsRotated=rotateTheta(coords2D,Theta)
#        # bringing the rotated coordinates into suited format (see above)
#        xyz3D=[[k[0],k[1],k[2]] for k in coordsRotated]
#
#        # Tritmp3: this will result into the traingulation of old and new plane
#        # Hence, the first 3 indices of each element triangulation need to be the same as ....
#        # the last 3 indices of each element of the former triangulation
#        # Former triangulation is given by all the Ntriangles last triangulationlist entries
#        Tritmp3 = Tri3[-Ntriangles:] # copying the triangutlation into Tritmp
#        Tritmp2 = [[kk[3], kk[4], kk[5]] for kk in Tritmp3] # setting the last 3 indices entries to the first three of the new triangulation
#        Tritmp3 = []
#
#        # Npoints is increasing with every additional plane
#        Npoints = len(coords3D)
#        # loop over all new elements/new triangulation
#        for kk in range(0,Ntriangles):
#            # append the temporary triangulation Tri3 element by element (same method as above)
#            # first add the first three indices (same as the last three from the prevoius triangulation slice)
#            Tritmp3.append(list(Tritmp2[kk]))
#            # then extend these 3 by 3 additional including the addition by the amount of points in the system (Npoints)
#            Tritmp3[kk].extend([kk+Npoints for kk in Tri2D[kk]])
#
#        #after finishing the new slices extend the final triangulation list
#        Tri3 += Tritmp3
#        # and the final coordinates list
#        coords3D += xyz3D
#    # after finishing for the whole 3D domain bring tringulation entries explicity into tineger format
#    # pyvtk will give an error outherwise
#    for k in range(0,len(Tri3)):
#        Tri3[k] = [int(kk) for kk in Tri3[k]]
#    return Tri3, coords3D
#
#
#def ParallelVideo(filename, i):
#    """
#    Generate a 3D visualization from 2D fields expanded azimuthally.
#    DEPRECATED - DOES NOT WORK ANYMORE
#
#    This function is used in parallel processes to export multiple frames of a video visualization in VTK format.
#    It is called within ExportSolutions.py and starts 50 processes (for 50 snapshots) or as many as possible.
#    One process serves for saving one scene.
#
#    Parameters
#    ----------
#    filename : str
#        Base name for loading input data and saving output files.
#    i : int
#        Index of the current snapshot for which the VTK file is created.
#    """
#    import pyvtk
#    # previously saved dictionary. Being loaded in each process to avoid conflicts or waiting times when sharing memory
#    # parallelization over scenes existing, meanes paralellization paramater is "i" in angle[i]
#    ExportDict = np.load(filename+'_Dict.npy').item()
#    # structure previously built in ExportSolution.py (structure for vtk export)
#    structure = ExportDict['structure']
#    #vtk point data init.
#    data=pyvtk.PointData()
#
#    #azimuthal wave number
#    m = ExportDict['m']
#
#    # amount of scene shots
#    n_Snaps=50
#    # dangle: increment in phase shift -> angle: including vector of all phase shifts occuring
#    dangle=2*np.pi/(n_Snaps-1)
#    dangle = np.arange(0,n_Snaps)*dangle
#    # for the 3D export also the phase in azimuthal direction is required <- 2pi/angularSteps
#    angularSteps = 100
#    # Init. of export dictionary
#    Export3D = {}
#    # theta is phase in azimuthal direction
#    theta = list(2*np.pi/angularSteps*np.arange(0,angularSteps+1))
#    # loop over all dictionary entries from imported dictionary, Howver: only entries with a ang-entry will be considered
#    for key in ExportDict:
#        if key[-3:] == 'ang': # if key has an angle entry continue
#            angle = ExportDict[key] # get phase field
#            Export3D[key[:-4]] = [] # init final dict entry for final field
#            # expand fields into azimuthal direction
#            # m             = azimuthal wave number
#            # theta[k]  = azimuthal angle of kth plane
#            # angle0    = unshifted phase angle (result from resolvent analyis)
#            # angle[i]  = temporal time shift with respect to the snapshot/process
#            # absol     = currently considered field's absolute values
#            # ang       = angle (previously computed)
#            # loop over azimuthally expanded planes
#            # each plane is the sum of phase field from 2D Felics plane
#            #       + phase due to azimuthal expansion
#            #       + phase due to temporal shift
#            # the final field results from projecting the field's absolute value onto the angle-field
#            for k in range(0,angularSteps+1):
#                angleNew = [m*theta[k]+angle0+dangle[i] for angle0 in ExportDict[key]]
#                Export3D[key[:-4]].extend( [ absol*np.cos(ang) for absol,ang in zip(ExportDict[key[:-4]+'_abs'],angleNew)])
#            data.append(pyvtk.Scalars(Export3D[key[:-4]], name=key[:-4]))
#
#    vtk = pyvtk.VtkData(structure, data)
#    vtk.tofile(filename+'Video_'+str(i)+'.vtk','binary')
#
#
#def smoothFieldWithKernel(field,FEMSpaces, iterations):
#    """
#    Smooth a scalar field on a finite element mesh using k-nearest neighbors. 
#    DEPRECATED. TO REVIVE ADD A WORKING TEST FOR IT IN TESTS
#
#    Parameters
#    ----------
#    field : dolfinx.fem.Function
#        The finite element field to smooth.
#    FEMSpaces : dict
#        Dictionary of FEM spaces (currently unused).
#    iterations : int
#        Number of smoothing iterations.
#
#    Returns
#    -------
#    dolfinx.fem.Function
#        The smoothed finite element field.
#    """
#
#    #from fenics import Function
#    from sklearn.neighbors import KNeighborsRegressor, RadiusNeighborsRegressor
#    smoothed_field = Function(field.function_space())
#    coords = field.function_space().tabulate_dof_coordinates()
#    neighbors = KNeighborsRegressor(n_neighbors=20, weights='uniform')
#    tmp = field.vector()[:]
#    for k in range(0,iterations):
#        neighbors.fit(coords,tmp)
#        tmp =  neighbors.predict(coords)
#
#    smoothed_field.vector()[:] = tmp
#
#    return smoothed_field
#
#
#def SutherlandLaw(param,MeanFlowDict):
#    """
#    Compute the molecular viscosity using Sutherland's law.
#
#    Parameters
#    ----------
#    param : object
#        Object containing viscosity and case parameters.
#    MeanFlowDict : dict
#        Dictionary of mean flow fields.
#
#    Returns
#    -------
#    np.ndarray
#        Array of kinematic viscosity values.
#    """
#    M=28.949
#    p=101300
#    R=8314.4598/M
#    mu0=param.Case.MolVisc
#    T0=291.15
#    C=120
#    if 'T' in param.Case.getMeanFlowFieldNames():
#        T = MeanFlowDict['T']
#    elif 'rho' in param.Case.getMeanFlowFieldNames():
#        T=p/R/MeanFlowDict['rho'].vector()[:]
#    else:
#        T = 293
#    nuMol=mu0*((T0+C)/(T+C)) * (T/T0)**1.5
#    return nuMol
#
#
#def executeParallelInterpolation( solutDir, nCubes):
#    """
#    Interpolate solution data onto finite element degrees of freedom for one cube in parallel.
#
#    Parameters
#    ----------
#    solutDir : str
#        Directory containing source data and saving destination pickle.
#    nCubes : int
#        Index of the current cube to interpolate.
#
#    Returns
#    -------
#    str
#        Path to the pickle file containing the interpolated values.
#    """
#    print (multiprocessing.current_process())
#    cubeTime = time.time()
#
#    # open the source dict and load the data
#    file = open(solutDir+'/SourceDict','rb')
#    data=pickle.load(file)
#    file.close()
#
#    # open the dofDict and load the data
#    file = open(solutDir+'/DofDict','rb')
#    DOFdata=pickle.load(file)
#    file.close()
#
#    print('min X' + str(min(data['pointsForCube'+str(nCubes)][:,0])))
#    print('max X' + str(max(data['pointsForCube'+str(nCubes)][:,0])))
#
#    print('min Y' + str(min(data['pointsForCube'+str(nCubes)][:,1])))
#    print('max Y' + str(max(data['pointsForCube'+str(nCubes)][:,1])))
#
#    print('min Z' + str(min(data['pointsForCube'+str(nCubes)][:,2])))
#    print('max Z' + str(max(data['pointsForCube'+str(nCubes)][:,2])))
#
#    # init destination Dict - each subprocess has its own destination dict (each subcube) - will be saved underneath
#    DestinationDict = {}
#    print(nCubes)
#    print(data['pointsForCube'+str(nCubes)].shape)
#    # this is the actual linear interpolation
#    # from (data['pointsForCube'+str(nCubes)], data['valsForCube'+str(nCubes)]) to DOFdata['dofCoordsOfCube'+str(nCubes)]
#    # the problem with linear interpolation is that if Felics points are exactly or slightly outside of origin grid
#    # there is no triangle or trapecoid which can be build around the Felics point --> NaN
#    # avoid NaN by fill_value=0 under the assumpotion that those boundary points are at a wall
#    DestinationDict['valsForCube'+str(nCubes)] = griddata(data['pointsForCube'+str(nCubes)],
#                                                          data['valsForCube'+str(nCubes)],
#                                                          DOFdata['dofCoordsOfCube'+str(nCubes)],
#                                                          method='linear', fill_value=0)
#
#    # under this name, the destination dict will be saved for each process/subcube and late loaded
#    pickleName = solutDir+'/pickeDumpForCube'+str(nCubes)
#    file = open(pickleName,'wb')
#    pickle.dump(DestinationDict,file)
#    file.close()
#
#    print('time for cube ' + str(time.time()-cubeTime))
#    return pickleName
#
#
#def parallelInterpolation(cubeBounds, dof_expanded, points, vals, tolerance, solutDir):
#    """
#    Perform parallel interpolation of volumetric data onto an expanded finite element grid.
#
#    Parameters
#    ----------
#    cubeBounds : np.ndarray
#        Boundaries of the cubes used for parallel interpolation.
#    dof_expanded : np.ndarray
#        Coordinates of the degrees of freedom in the expanded grid.
#    points : np.ndarray
#        Coordinates of the original data points.
#    vals : np.ndarray
#        Values associated with the original data points.
#    tolerance : float
#        Tolerance to expand source cubes to ensure coverage.
#    solutDir : str
#        Directory path for storing temporary pickles.
#
#    Returns
#    -------
#    np.ndarray
#        Interpolated values on the expanded grid.
#    """
#    import time
#    import os
#
#    dofDict = {} # is passed to the pool-processes later: dof coordinates of azimuthally expanded Felics grid
#    SourceDict = {} # is passed to the pool-process later: origin data of origin sub cubes
#    for cubeCnt in range(0,cubeBounds.shape[0]): # loop over all sub cubes built previously
#        # find indices of expanded Felics grid within sub cubes of origin domain
#        # in x,y,z each, then, intersection of them are the coordinates indices of inside the subcubes
#        Xidx = np.argwhere((dof_expanded[:,0]>=cubeBounds[cubeCnt,0]) & (dof_expanded[:,0]<=cubeBounds[cubeCnt,1])).flatten()
#        Yidx = np.argwhere((dof_expanded[:,1]>=cubeBounds[cubeCnt,2]) & (dof_expanded[:,1]<=cubeBounds[cubeCnt,3])).flatten()
#        Zidx = np.argwhere((dof_expanded[:,2]>=cubeBounds[cubeCnt,4]) & (dof_expanded[:,2]<=cubeBounds[cubeCnt,5])).flatten()
#        IDX = np.intersect1d(np.intersect1d(Xidx, Yidx),Zidx)
#        # dof Dict- incorporates the indices for each subcube of the expanded Felics grid
#        dofDict['dofIDXofCube'+str(cubeCnt)] = IDX
#        # ....and their coordinates
#        dofDict['dofCoordsOfCube'+str(cubeCnt)] = dof_expanded[IDX,:]
#
#        # get the bounds of the currently considered origin subcube and apply tolerance band such that
#        # the origin subcubes become slightly larger than destination subcubes to avoid missed coordinates
#        sourceCube = [cubeBounds[cubeCnt,0]-tolerance, cubeBounds[cubeCnt,1]+tolerance,\
#                        cubeBounds[cubeCnt,2]-tolerance, cubeBounds[cubeCnt,3]+tolerance,\
#                        cubeBounds[cubeCnt,4]-tolerance, cubeBounds[cubeCnt,5]+tolerance ]
#
#        # find all coordinates of the origin data within the subcube boundaries inclusive the tolerance
#        Xidx = np.argwhere((points[:,0]>=sourceCube[0]) & (points[:,0]<=sourceCube[1])).flatten()
#        Yidx = np.argwhere((points[:,1]>=sourceCube[2]) & (points[:,1]<=sourceCube[3])).flatten()
#        Zidx = np.argwhere((points[:,2]>=sourceCube[4]) & (points[:,2]<=sourceCube[5])).flatten()
#        IDX = np.intersect1d(np.intersect1d(Xidx, Yidx),Zidx)
#        # save coordiantes and values of the origin sub cubes in source dict for each cube (cubeCNT)
#        SourceDict['pointsForCube'+str(cubeCnt)] = points[IDX,:]
#        SourceDict['valsForCube'+str(cubeCnt)] = vals[IDX,:]
#
#    # method to pass data to subprocesses. after a certain data size the multiprocessing library fails to pass them as args
#    # time measure shows that this is very cheap
#    # dicts are save in the same directory as Felics output
#    pickleDumpTime = time.time()
#    file = open(solutDir+'/SourceDict','wb')
#    pickle.dump(SourceDict,file)
#    file.close()
#    print('time for pickleDump:' +str(time.time()-pickleDumpTime))
#
#    pickleDumpTime = time.time()
#    file = open(solutDir+'/DofDict','wb')
#    pickle.dump(dofDict,file)
#    file.close()
#    print('time for pickleDump:' +str(time.time()-pickleDumpTime))
#
#    print('Starting parallel interpolation....')
#    # a pool of 12 processes HARD CODED CURRENTLY
#    pool=multiprocessing.Pool(12)
#    # execute the executeParallelInterpolation-function on each proces
#    func= partial(executeParallelInterpolation, solutDir)
#    # the variation for each proces is a simple coutner from 0 to length(subCubes)
#    # the function executeParallelInterpolation loads the above pickled data
#    # .... interpolates each subcubes
#    # .... and pickles the data of each sub-cube again
#    # .... the file name of the interpolated and pickled destination dat is returned to pickleNames
#    pickleNames=pool.map(func, range(0,cubeBounds.shape[0]))
#
#    print('Done with parallel processes')
#
#
#    # the destination subcubes, i.e. the pickled destination data from executeParallelInterpolation , are assembled to global cube
#    print('starting to assemble subcubes')
#    assembleTime = time.time()
#    # the shape of the destination data is given by the amount of dof-coordinates of expanded Felcis grid
#    # .... and the amount of quanties (e.g. in P2 interpoaltion on inerpolates U,V,W simultanously)
#    newVals = np.zeros((dof_expanded.shape[0],vals.shape[1]))
#    for cubeCnt in range(0,cubeBounds.shape[0]):
#        pickeLoadTime = time.time()
#        # open pickle file
#        file = open(pickleNames[cubeCnt],'rb')
#        # get data instisde opened pickle file
#        data = pickle.load(file)
#        print('time to pickle load: '+str(time.time()-pickeLoadTime))
#        # here is the reason why we save the indices of the Felics subcube coordinates relative to the global dof-list
#        # we need to point to the correct place in order to assemble the subcubes
#        # dofDict['dofIDXofCube'+str(cubeCnt)] ahs the indices of each subcube in the global dof-list
#        newVals[dofDict['dofIDXofCube'+str(cubeCnt)], :] = data['valsForCube'+str(cubeCnt)]
#        # remove the file which was save by the pool-process previously
#        os.remove(pickleNames[cubeCnt])
#    print('returning interpolated vals...assembling took: '+str(assembleTime))
#    return newVals
#
#def loadCSV(path):
#    """
#    Load a CSV file into a dictionary with numpy arrays for each column.
#
#    Parameters
#    ----------
#    path : str
#        Path to the CSV file.
#
#    Returns
#    -------
#    dict
#        Dictionary containing column headers as keys and their values as arrays.
#    """
#    import numpy as np
#    import csv
#    with open(path, newline='') as f:
#        reader = csv.reader(f)
#        row_count = sum(1 for row in reader)
#    with open(path, newline='') as f:
#        reader = csv.reader(f)
#        index=0
#        for row in reader:
#            if index==0:
#                TableDict={}
#                header=row
#                for entry in header:
#                    TableDict[entry]=np.zeros(row_count-1)
#            else:
#                for i in range(0,len(row)):
#                    TableDict[header[i]][index-1]=row[i]
#            index+=1
#
#    return TableDict
