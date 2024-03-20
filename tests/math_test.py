# This is for math test

import numpy as np

# Test vector projection onto a plane
# normal vector of the plane: n = [0, 0, 1]
# vector to be projected: v = [1, 1, 1]
# projection of v onto the plane: v - (v.n)*n = [1, 1, 0]

n = np.array([0, 0, 1])
v = np.array([1, 1, 1])
v_proj = v - np.dot(v, n)*n

print("Projection of v onto the plane = \n{:}\n".format(v_proj))
