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
import numpy as np
import sys

# call function in terminal via: "python attila2felics.py path/to/mesh.inp"

# functions
def find_line_number(lines, search_string):  # finds and returns number of a line that starts with given string
    for current_line_number, line in enumerate(lines, start=1):
        if line.startswith(search_string):
            return current_line_number
    return None  # if no line starts with the given string

# main code
if len(sys.argv) == 1:
    icemFile = input("Please provide path of .inp mesh file as argument!")
else:
    icemFile = sys.argv[1]

# Read the entire file into memory
with open(icemFile, 'r') as file1:
    lines = file1.readlines()

# Search for the number of dimensions, elements, vertices, and sides in a single loop
ndim = nelem = nvert = nsides = None
for i, line in enumerate(lines):
    if line.startswith('  ndim'):
        ndim = int(line[18:])
    elif line.startswith('  ncells'):
        nelem = int(line[18:])
    elif line.startswith('  nnodes '):
        nvert = int(line[18:])
    elif line.startswith('  nsides '):
        nsides = int(line[18:])

if None in (ndim, nelem, nvert, nsides):
    print("Error: Missing mesh metadata")
    sys.exit(1)

print(f"Number of mesh dimensions: {ndim}")
print(f"Number of mesh elements: {nelem}")
print(f"Number of mesh vertices: {nvert}")
print(f"Number of sides: {nsides}")

# Reading vertices
print('Reading vertices...')
vert_line_number = find_line_number(lines, 'nodes')  # find start of node listing
vert = np.zeros((nvert, ndim))

for i in range(vert_line_number, vert_line_number + nvert):
    line_entries = lines[i].split()  # Split line into max 4 parts
    vert[i - vert_line_number, :ndim] = np.array(line_entries[1:4], dtype=float) 

vert /= 1000    # Convert to meters

# Finding number of boundaries
nmark_line_number = find_line_number(lines, 'end_side_flags')
nmark = nmark_line_number - find_line_number(lines, 'side_flags') - 2
print(f"Number of boundaries: {nmark}")

# Read BCs
print('Reading boundary conditions...')
BCs = np.zeros(nmark, dtype=int)
BC_names = []
BC_line_number = find_line_number(lines, 'side_flags')
for i in range(BC_line_number + 1, nmark_line_number - 1):
    line_entries_bc = lines[i].split()
    BCs[i - BC_line_number - 1] = int(line_entries_bc[0])
    BC_names.append(line_entries_bc[1])
print(BCs)

# Reading cells
print('Reading cells...')
cell_line_number = find_line_number(lines, 'cells')     # find start of cell listing
cells = np.zeros((nelem, ndim+1), dtype=int)
for i in range(cell_line_number, cell_line_number + nelem):
    line_entries = lines[i].split()                     # Split line into max 6 parts
    cells[i - cell_line_number, :ndim+1] = np.array(line_entries[2:6], dtype=int)

# Reading sides
print('Reading sides...')
side_line_number = find_line_number(lines, 'sides')     # find start of side listing
sides = np.zeros((nsides, ndim+1), dtype=int)
for i in range(side_line_number, side_line_number + nsides):
    line_entries = lines[i].split()
    sides[i - side_line_number, :ndim+1] = np.array(line_entries[2:6], dtype=int)
    for j in range(nmark):
         if sides[i - side_line_number,3] == BCs[j]:
              sides[i - side_line_number,3] = j + 1

# Sorting sides by boundary condition index
print('Sorting sides...')
counter_sort = 0
sides_sort = np.zeros((nsides, ndim+1))
for j in range(1, nmark+1):
     for i in range(0, nsides):
          if sides[i,3] == j:
               sides_sort[counter_sort,:] = sides[i,:]
               counter_sort = counter_sort + 1       
sides = sides_sort.astype(int)

# Writing the new mesh file
exportfilename = icemFile[0:-4] + '.msh'
print(f'Exporting mesh to {exportfilename}')
with open(exportfilename, 'w') as file2:
    file2.write('$MeshFormat\n')
    file2.write('2.2 0 8\n')
    file2.write('$EndMeshFormat\n')
    file2.write('$PhysicalNames\n')
    file2.write(f'{nmark}\n')
    
    # Write boundary conditions with numbers and names
    for i in range(nmark):
        file2.write(f'1 {i + 1} "{BC_names[i]}"\n')
    
    file2.write('$EndPhysicalNames\n')
    file2.write('$Nodes\n')
    file2.write(f'{nvert}\n')
    for i in range(nvert):
        file2.write(f'{i+1} {vert[i, 0]} {vert[i, 1]} {vert[i, 2]}\n')
    
    file2.write('$EndNodes\n')
    file2.write('$Elements\n')
    file2.write(f'{nelem + nsides}\n')

    # Writing side elements
    elem_type = {2: '1', 3: '2'}[ndim]  # Dictionary for element types
    for i in range(nsides):
        file2.write(f'{i+1} {elem_type} 2 {sides[i, 3]} {sides[i, 3]} {sides[i, 0]} {sides[i, 1]} {sides[i, 2]}\n')
    
    # Writing cell elements
    elem_type = {2: '2', 3: '4'}[ndim]  # Adjusted element type for cells
    for i in range(nelem):
        file2.write(f'{i + nsides + 1} {elem_type} 2 {nmark + 1} {nmark + 1} {cells[i, 0]} {cells[i, 1]} {cells[i, 2]} {cells[i, 3]}\n')
    
    file2.write('$EndElements\n')
