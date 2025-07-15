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


import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('TkAgg')

x_coords, y_coords = read_csv_file("./out/spectrum.csv")
# Plotting
fig = plt.figure(figsize=(9, 6))
plt.rcParams.update({'font.size': 24})
plt.scatter(x_coords, y_coords, c='blue', marker='o')
plt.title('Scatter Plot of Eigenvalues')
plt.xlabel('Im(Eigenvalue)')
plt.ylabel('Re(Eigenvalue)')
plt.grid(True)
plt.show()
