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
'''This code transforms su2 mesh files produced for example by the centaur software to .msh files, to be used with FELiCS'''
'''Author: T.L.Kaiser'''
import numpy as np
import sys 

if len(sys.argv)==1:
	su2File=input("Please provide path of su2 mesh file as argument!")
else:
	su2File= sys.argv[1]
# Using readlines() 
file1 = open(su2File, 'r') 

line=file1.readline()
if line[0:6] == 'NDIME=':
	ndim=int(line[7:])
	print('Number of mesh dimensions: '+str(ndim))
else:
	print('Error! Number of mesh dimensions not found')

line=file1.readline()
if line[0:6] == 'NELEM=':
	nelem=int(line[7:])
	print('Number of mesh elements: '+str(nelem))
else:
	print('Error! Number of elements not found')
print('Reading triangulation...')
tri=np.zeros((nelem,3))
for i in range(0,nelem):
	line=file1.readline()[2:]
	tri[i,:]=np.array(line.split('  ')[1:4])

line=file1.readline()
if line[0:6] == 'NPOIN=':
	nvert=int(line[7:])
	print('Number of mesh vertices: '+str(nvert))
else:
	print('Error! Number of vertices not found')
print('Reading vertices...')
vert=np.zeros((nvert,2))
for i in range(0,nvert):
	line=file1.readline()
	vert[i,:]=np.array(line.split('  ')[0:2])


line=file1.readline()
if line[0:6] == 'NMARK=':
	nmark=int(line[7:])
	print('Number of boundaries: '+str(nmark))
else:
	print('Error! Number of boundaries not found')
boundaries=list()
print('Reading boundaries...')
nelem_tot=nelem+0
for mark_i in range(0,nmark):
	boundary_dict={}
	line=file1.readline()
	if line[0:11] == 'MARKER_TAG=':
		boundary_dict['name']=line[12:].strip()
	else:
		print('Error! Boundary name not found')
	line=file1.readline()
	if line[0:13] == 'MARKER_ELEMS=':
		boundary_dict['nelem']=int(line[14:])
		nelem_tot+=boundary_dict['nelem']
	else:
		print('Error! number of Boundary elements not found')
	boundary_dict['elements']=np.zeros((boundary_dict['nelem'],2))
	print('\tReading '+str(boundary_dict['nelem'])+' elements of '+boundary_dict['name']+'...')
	boundary_dict['elements']=np.zeros((nelem,2))
	for i in range(0,boundary_dict['nelem']):
		line=file1.readline()[2:]
		boundary_dict['elements'][i,:]=np.array(line.split('  ')[1:3])
	boundaries.append(boundary_dict)

exportfilename=su2File[0:-8]+'.msh'
print('Exporting mesh to '+exportfilename)
file2=open(exportfilename, 'w')
file2.write('$MeshFormat\n')
file2.write('2.2 0 8\n')
file2.write('$EndMeshFormat\n')
file2.write('$PhysicalNames\n')
file2.write(str(len(boundaries)+1)+'\n')
for i in range(0,len(boundaries)):
	file2.write('1 '+str(i+1)+' \"'+boundaries[i]['name']+'\"'+'\n')
file2.write('2 '+str(i+2)+' \"all\"'+'\n')
file2.write('$EndPhysicalNames\n')
file2.write('$Nodes\n')
file2.write(str(nvert)+'\n')
for i in range(0,nvert):
	file2.write(str(i+1) +' '+ str(vert[i,0])+' '+str(vert[i,1])+ ' 0\n')
file2.write('$EndNodes\n')
file2.write('$Elements\n')
file2.write(str(nelem_tot)+'\n')
i_inter=0
for boundary_index, boundary in enumerate(boundaries):
	print(boundary['nelem'])
	for i in range(0,boundary['nelem']):
		file2.write(str(i+1+i_inter)+ ' 1 2 '+str(boundary_index+1)+ ' ' +str(boundary_index+1)+' '+str(int(boundary['elements'][i,0])+1)+ ' ' +str(int(boundary['elements'][i,1])+1)+'\n')
	i_inter+=i
for i in range(0,nelem):
	file2.write(str(i+1+i_inter)+ ' 2 2 '+str(nmark+1)+ ' ' +str(nmark+1)+' '+str(int(tri[i,0])+1)+ ' ' +str(int(tri[i,1])+1)+' '+str(int(tri[i,2])+1)+' \n')
file2.write('$EndElements\n')

