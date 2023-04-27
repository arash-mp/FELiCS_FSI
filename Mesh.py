#/* Copyright (C) 2019 FLOW group TU Berlin - All Rights Reserved
# * You may NOT use, distribute or modify this code without explicit
# * permission of the copyright owner, the FLOW group at TU Berlin!
# * However, permissions to use and modify the code are generally 
# * granted when asked for.
# * To ask for permission please contact t.kaiser@tu-berlin.de.
# */
import numpy as np
from fenics import Mesh
import os
#import meshio
'''
# **********************************************************************
# * This file provides a routine for reading the mesh
# * 
# * This file was created by Thomas L. Kaiser. Significant contributions 
# * were made by
# * - 
# *   
# *   
# ********************
'''

def ReadMesh(path):
	''' Function reading mesh from path "path"
	\t Possible mesh types: 
	\t \t -xml'''
	if path.split('.')[-1]=='xml':
		## Read mesh from xml file. 
		print ('Reading mesh in XML format from "'+path+'"...')
		mesh=Mesh(path)
	elif path.split('.')[-1]=='msh':
#		print ('Convert mesh from .msh to .xml using dolfin-convert gmshnetzname.msh dolfinmeshname.xml in "'+path+'"...')
		print('Converting mesh from msh(gmsh) format to .xml')
		os.system('dolfin-convert '+path+' '+path[:-4]+'.xml')
		mesh = Mesh(path[:-4]+'.xml')
		
        
        
        
	elif path.split('.')=='dat':
		print ('Reading mesh in .dat format from "'+path+'"...')
		#'mesh=Mesh(path)
	return mesh
