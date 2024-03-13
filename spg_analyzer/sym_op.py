from calendar import c
import numpy as np
# symmetry operation functions

from data import TOLERANCE

# Inversion



# Proper Rotation

# Cn

def Cn_rot(coord, n, xyz):
    '''
    Rotate the coordinates by Cn.

    Args:
        coord: coordinates of the molecule
        n: order of rotation, 2, 3, 4, 5, 6, ...
        xyz: axis of rotation, 'x', 'y', or 'z'

    Rotational matrix:
    z-axis:
    [[cos(theta), -sin(theta), 0],
    [sin(theta),   cos(theta), 0],
    [     0,          0,       1]]

    y-axis:
    [[cos(theta), 0, sin(theta)],
    [     0,      1,     0     ],
    [-sin(theta), 0, cos(theta)]]

    x-axis:
    [[1,    0,            0    ],
    [0, cos(theta), -sin(theta)],
    [0, sin(theta), cos(theta)]]

    where theta = 2*pi/n

    '''
    rot_mat = np.zeros((3,3))
    theta = 2*np.pi/n
    rot_mat[0][0] = np.cos(theta)
    rot_mat[0][1] = -np.sin(theta)
    rot_mat[1][0] = np.sin(theta)
    rot_mat[1][1] = np.cos(theta)
    rot_mat[2][2] = 1

    if xyz == 'z':
        rot_mat = rot_mat
        # print("rotation matrix z", rot_mat) 
    elif xyz == 'y':
        rot_mat = rot_mat.T
        rot_mat[[1,2],:] = rot_mat[[2,1],:]
        rot_mat[:,[1,2]] = rot_mat[:,[2,1]]
        # print("rotation matrix y", rot_mat) 
    elif xyz == 'x':
        rot_mat[[0,1,2],:] = rot_mat[[2,0,1],:]
        rot_mat[:,[0,1,2]] = rot_mat[:,[2,0,1]]
        # print("rotation matrix x", rot_mat)
    else:
        raise ValueError('Invalid axis of rotation.')

    return np.dot(rot_mat, coord.T).T 


# Improper Rotation

# Reflection

def reflection(coord, xyz):
    '''
    Reflect the coordinates.

    Args:
        coord: coordinates of the molecule
        xyz: plane of reflection, 'xy', 'yz', or 'xz'

    Reflection matrix:
    xy-plane:
    [[1, 0, 0],
    [0, 1, 0],
    [0, 0, -1]]

    yz-plane:
    [[-1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]]

    xz-plane:
    [[1, 0, 0],
    [0, -1, 0],
    [0, 0, 1]]

    '''

    ref_mat = np.zeros((3,3))
    if xyz == 'xy':
        ref_mat[0][0] = 1
        ref_mat[1][1] = 1
        ref_mat[2][2] = -1
    elif xyz == 'yz':
        ref_mat[0][0] = -1
        ref_mat[1][1] = 1
        ref_mat[2][2] = 1
    elif xyz == 'xz':
        ref_mat[0][0] = 1
        ref_mat[1][1] = -1
        ref_mat[2][2] = 1
    else:
        raise ValueError('Invalid plane of reflection.')
    
    return np.dot(ref_mat, coord.T).T


