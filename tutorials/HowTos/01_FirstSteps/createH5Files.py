import numpy as np
import h5py

xx = np.linspace(-0.1, 1.1, 101)
yy = np.linspace(-0.1, 1.1, 101)
X, Y = np.meshgrid(xx, yy, indexing='ij')
x = (X.flatten())
y = (Y.flatten())


function_sine = np.sin(2.*np.pi*x) * np.sin(2.*np.pi*y)

file = "function_sine.h5"

with h5py.File(file, 'w') as f:
    f.create_dataset('x',   data=x)
    f.create_dataset('y',   data=y)
    f.create_dataset('function_sine', data = function_sine)
