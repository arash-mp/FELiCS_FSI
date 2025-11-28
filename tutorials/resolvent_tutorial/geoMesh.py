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
# Import necessary libraries
import numpy as np

# Define the path and name for the output .geo file
file_name = 'FeliCS_mesh.geo'

# Load the upper boundary points from LES data and normalize by diameter
upB2D = np.load('upB2D.npy')

# Calculate maximum dimensions from the boundary points
MAX_x = np.max(upB2D[:, 0])
MAX_r = np.max(upB2D[:, 1])

# Mesh generation parameters
CellsFineness = 0.5  # Scaling factor for number of cells
mesh_size = 0.01  # Base mesh size
Nxfactor = 0.5  # Factor to adjust the number of cells in the x-direction

# Create and write the .geo file
with open(file_name, 'w') as geo_file:
    # Define the base mesh size
    geo_file.write(f'sm_top = {mesh_size};\n')

    # Write the wall points (upper boundary)
    for i, (x, y) in enumerate(upB2D):
        geo_file.write(f'Point({i+1}) = {{{x}, {y}, 0, sm_top}};\n')

    # Define key points along the bottom boundary
    x_mid = 100  # Midpoint index
    x_left = 80  # Left segment boundary
    x_right = 150  # Right segment boundary

    # Define bottom points
    bottom_points = [
        (len(upB2D) + 1, upB2D[0, 0], 0),  # Leftmost point
        (len(upB2D) + 2, upB2D[x_left, 0], 0),  # Left segment point
        (len(upB2D) + 3, upB2D[x_mid, 0], 0),  # Midpoint
        (len(upB2D) + 4, upB2D[x_right, 0], 0),  # Right segment point
        (len(upB2D) + 5, MAX_x, 0)  # Rightmost point
    ]

    # Write bottom points to the file
    for p_id, x, y in bottom_points:
        geo_file.write(f'Point({p_id}) = {{{x}, {y}, 0, sm_top}};\n')

    # Define the spline segments for the upper wall
    for segment, (start, end) in enumerate(zip([0, x_left, x_mid, x_right],
                                             [x_left, x_mid, x_right, len(upB2D)-1])):
        point_indices = ', '.join(str(i+1) for i in range(start, end + 1))
        geo_file.write(f'Spline({segment+1}) = {{{point_indices}}};\n')

    # Define vertical and horizontal lines to complete the domain
    segment = 4  # Current segment counter
    geo_file.write(f'Line({segment+1}) = {{{len(upB2D)+5}, {len(upB2D)}}};\n')  # Right vertical line
    geo_file.write(f'Line({segment+2}) = {{{len(upB2D)+1}, 1}};\n')  # Left vertical line

    # Define bottom horizontal lines
    geo_file.write(f'Line({segment+3}) = {{{len(upB2D)+5}, {len(upB2D)+4}}};\n')
    geo_file.write(f'Line({segment+4}) = {{{len(upB2D)+4}, {len(upB2D)+3}}};\n')
    geo_file.write(f'Line({segment+5}) = {{{len(upB2D)+3}, {len(upB2D)+2}}};\n')
    geo_file.write(f'Line({segment+6}) = {{{len(upB2D)+2}, {len(upB2D)+1}}};\n')

    # Define vertical connectors between top and bottom
    geo_file.write(f'Line({segment+7}) = {{{len(upB2D)+2}, {x_left+1}}};\n')
    geo_file.write(f'Line({segment+8}) = {{{len(upB2D)+3}, {x_mid+1}}};\n')
    geo_file.write(f'Line({segment+9}) = {{{len(upB2D)+4}, {x_right+1}}};\n')

    # Define mesh refinement parameters for different regions
    num_cells_v, progress_v = 50, 1.03  # Vertical mesh parameters
    progress_vS4 = 1.1  # Special progression for outlet region

    # Region-specific parameters
    progress_S1, progress_S2, progress_S3, progress_S4 = 1+0.014/Nxfactor, 1+0.04/Nxfactor, 1, 1/(1 + 0.018/Nxfactor)
    num_cells_S1, num_cells_S2, num_cells_S3, num_cells_S4 = int(Nxfactor*60), int(Nxfactor*40), int(Nxfactor*300), int(Nxfactor*200)

    # Apply transfinite mesh to each region
    # Region 1
    geo_file.write(f'Transfinite Curve{{{segment+2}}} = {num_cells_v} Using Progression {1/progress_v};\n')  # Left line
    geo_file.write(f'Transfinite Curve{{{segment+7}}} = {num_cells_v} Using Progression {1/progress_v};\n')  # Right vertical
    geo_file.write(f'Transfinite Curve{{1}} = {num_cells_S1} Using Progression {1/progress_S1};\n')  # Top
    geo_file.write(f'Transfinite Curve{{{segment+6}}} = {num_cells_S1} Using Progression {progress_S1};\n')  # Bottom

    # Region 2
    geo_file.write(f'Transfinite Curve{{{segment+5}}} = {num_cells_S2} Using Progression {progress_S2};\n')
    geo_file.write(f'Transfinite Curve{{2}} = {num_cells_S2} Using Progression {1/progress_S2};\n')  # Top
    geo_file.write(f'Transfinite Curve{{{segment+8}}} = {num_cells_v} Using Progression {1/progress_v};\n')

    # Region 3
    geo_file.write(f'Transfinite Curve{{{segment+9}}} = {num_cells_v} Using Progression {1/progress_v};\n')
    geo_file.write(f'Transfinite Curve{{3}} = {num_cells_S3} Using Progression {1/progress_S3};\n')  # Top
    geo_file.write(f'Transfinite Curve{{{segment+4}}} = {num_cells_S3} Using Progression {progress_S3};\n')

    # Region 4
    geo_file.write(f'Transfinite Curve{{{segment+3}}} = {num_cells_S4} Using Progression {progress_S4};\n')
    geo_file.write(f'Transfinite Curve{{4}} = {num_cells_S4} Using Progression {1/progress_S4};\n')  # Top
    geo_file.write(f'Transfinite Curve{{{segment+1}}} = {num_cells_v} Using Progression {1/progress_vS4};\n')

    # Define surface loops and apply transfinite meshing
    geo_file.write('Line Loop(6) = {1, -11, 10, 6};\n')
    geo_file.write('Plane Surface(7) = {6};\n')
    geo_file.write('Line Loop(8) = {2, -12, 9, 11};\n')
    geo_file.write('Plane Surface(9) = {8};\n')
    geo_file.write('Line Loop(10) = {3, -13, 8, 12};\n')
    geo_file.write('Plane Surface(11) = {10};\n')
    geo_file.write('Line Loop(12) = {4, -5, 7, 13};\n')
    geo_file.write('Plane Surface(13) = {12};\n')
    geo_file.write('Transfinite Surface {7, 9, 11, 13};\n')

    # Define physical groups for boundary conditions
    geo_file.write('Physical Surface("domain") = {7, 9, 11, 13};\n')
    geo_file.write('Physical Curve("centerline") = {7, 8, 9, 10};\n')
    geo_file.write('Physical Curve("inlet") = {6};\n')
    geo_file.write('Physical Curve("outlet") = {5};\n')
    geo_file.write('Physical Curve("wall") = {1, 2, 3, 4};\n')

print(f'Gmsh .geo file "{file_name}" has been successfully created.')
