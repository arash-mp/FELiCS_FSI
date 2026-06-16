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
"""
    Scatter plot from a CSV file. 
    Author      : [Anant Talasikar]
    Date        : [25.06.2025]

    Description:
    The function <read_csv_file>, reads a CSV file with headers 'omega_direct_r' and 'omega_direct_i' and 
    returns them as x_coords and y_coords.

    Parameters:
    file_path (str): The path to the CSV file.
    """

import csv

def read_csv_file(file_path):
    x_coords = []
    y_coords = []
    # print
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            headers = next(reader)  # Skip the header row
            
            for row in reader:
                if len(row) >= 2:
                    try:
                        x = float(row[0])
                        y = float(row[1])
                        x_coords.append(x)
                        y_coords.append(y)
                    except ValueError:
                        print(f"Skipping invalid row: {row}")
        return x_coords, y_coords
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

x_coords, y_coords = read_csv_file("./output_dir/spectrum.csv")
# Plotting
fig = plt.figure(figsize=(9, 6))
plt.rcParams.update({'font.size': 24})
plt.scatter(x_coords, y_coords, c='blue', marker='o')
plt.title('Scatter Plot of Eigenvalues')
plt.xlabel('Re(Eigenvalue)')
plt.ylabel('Im(Eigenvalue)')
plt.grid(True)
plt.savefig('eigenvalues_scatter_plot.png')
plt.close()
