# Improper Rotation
import numpy as np

from utils import sort_coords, thre_cut

def Sn_rot(coord, n, xyz):
    '''
    Rotate and reflects the coordinates by Sn.

    Args:
        coord: coordinates of the molecule
        n: order of rotation, 2, 3, 4, 5, 6, ...
        xyz: axis of rotation, 'x', 'y', or 'z'

    Rotational matrix:
    z-axis:
    [[cos(theta), -sin(theta), 0],  # rotation matrix
    [sin(theta),   cos(theta), 0],
    [     0,          0,       1]]

    y-axis:
    [[cos(theta), 0, sin(theta)],
    [     0,      1,     0     ],
    [-sin(theta), 0, cos(theta)]]  # rotation matrix

    x-axis:
    [[1,    0,            0    ],
    [0, cos(theta), -sin(theta)],
    [0, sin(theta), cos(theta]]]  # rotation matrix

    where theta = 2*pi/n

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
    
    rot_mat = np.zeros((3,3))
    ref_mat = np.zeros((3,3))
    theta = 2*np.pi/n
    rot_mat[0][0] = np.cos(theta)
    rot_mat[0][1] = -np.sin(theta)
    rot_mat[1][0] = np.sin(theta)
    rot_mat[1][1] = np.cos(theta)
    rot_mat[2][2] = 1

    if xyz == 'z':
        rot_mat = rot_mat
        ref_mat[0][0] = 1
        ref_mat[1][1] = 1
        ref_mat[2][2] = -1
        # print("rotation matrix z", rot_mat) 
    elif xyz == 'y':
        rot_mat = rot_mat.T
        rot_mat[[1,2],:] = rot_mat[[2,1],:]
        rot_mat[:,[1,2]] = rot_mat[:,[2,1]]
        ref_mat[0][0] = 1
        ref_mat[1][1] = -1
        ref_mat[2][2] = 1
        # print("rotation matrix y", rot_mat) 
    elif xyz == 'x':
        rot_mat[[0,1,2],:] = rot_mat[[2,0,1],:]
        rot_mat[:,[0,1,2]] = rot_mat[:,[2,0,1]]
        ref_mat[0][0] = 1
        ref_mat[1][1] = -1
        ref_mat[2][2] = 1
        # print("rotation matrix x", rot_mat)
    else:
        raise ValueError('Invalid axis of rotation.')
    irot_mat=np.dot(rot_mat, ref_mat)
    return np.dot(irot_mat, coord.T).T
    
def check_Sn(coord):
    '''
    Check the symmetry of the molecule using Sn (from S8 - S3), excluding S2 and S1.

    Args:
        coord: coordinates of the molecule

    Returns:
        order of rotation, n, if the molecule has Sn symmetry, None otherwise.
    '''

    for i in range(3, 9).__reversed__():
        irot_z = sort_coords(thre_cut(Sn_rot(coord, i, 'z'), 4))
        irot_y = sort_coords(thre_cut(Sn_rot(coord, i, 'y'), 4))
        irot_x = sort_coords(thre_cut(Sn_rot(coord, i, 'x'), 4))
        coord = sort_coords(thre_cut(coord, 4))
        # print(f"rot_z: {rot_z}")
        # print(f"rot_y: {rot_y}")
        # print(f"rot_x: {rot_x}")
        # print(f"coord: {coord}")

        if np.array_equal(irot_z, coord):
            print(f"S{i} improper rotation around z-axis is a symmetry operation for the molecule.")
            return i
        elif np.array_equal(rot_y, coord):
            print(f"S{i} improper rotation around y-axis is a symmetry operation for the molecule.")
            return i
        elif np.array_equal(rot_x, coord):
            print(f"S{i} improper rotation around x-axis is a symmetry operation for the molecule.")
            return i
        else:
            continue
    
    print("No Sn (3-8) rotation is a symmetry operation for the molecule.")
    return None