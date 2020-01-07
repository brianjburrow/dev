import numpy as np

from prismatic_beam import prismatic_element
n_Beams = 4 

orientations = np.array([90., 45., -45., -90.])

orientations *= np.pi / 180.0

A = 1
I = 1
L = 1
E = 1

elements = [prismatic_element(A, I, L, E, orientation) for orientation in orientations] 

for element in elements:
    print(element.rotation_matrix.T)

