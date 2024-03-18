import numpy as np
# symmetry operation functions

from utils import delta, sort_coords, thre_cut
from data import extend_axis

# Proper Rotation

# Cn

def Cn_rot(coord, n, xyz):
    '''
    Rotate the coordinates by Cn along x, y, or z axis.

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

    return np.matmul(rot_mat, coord.T).T 


def Cn_rot_axis(coord, n, axis):
    '''
    Rotate the coordinates by Cn along a given axis unit vector.

    Args:
        coord: coordinates of the molecule
        n: order of rotation, 2, 3, 4, 5, 6, ...
        axis: unit vector of axis of rotation, [ux, uy, uz], where ux^2 + uy^2 + uz^2 = 1

    Rotational matrix:
    [[cos(theta) + ux^2*(1-cos(theta)),    ux*uy*(1-cos(theta)) - uz*sin(theta), ux*uz*(1-cos(theta)) + uy*sin(theta)],
    [uy*ux*(1-cos(theta)) + uz*sin(theta), cos(theta) + uy^2*(1-cos(theta)),     uy*uz*(1-cos(theta)) - ux*sin(theta)],
    [uz*ux*(1-cos(theta)) - uy*sin(theta), uz*uy*(1-cos(theta)) + ux*sin(theta), cos(theta) + uz^2*(1-cos(theta))]]

    where theta = 2*pi/n

    '''
    rot_mat = np.zeros((3,3))
    theta = 2*np.pi/n
    ux, uy, uz = axis
    rot_mat[0][0] = np.cos(theta) + ux**2*(1-np.cos(theta))
    rot_mat[0][1] = ux*uy*(1-np.cos(theta)) - uz*np.sin(theta)
    rot_mat[0][2] = ux*uz*(1-np.cos(theta)) + uy*np.sin(theta)
    rot_mat[1][0] = uy*ux*(1-np.cos(theta)) + uz*np.sin(theta)
    rot_mat[1][1] = np.cos(theta) + uy**2*(1-np.cos(theta))
    rot_mat[1][2] = uy*uz*(1-np.cos(theta)) - ux*np.sin(theta)
    rot_mat[2][0] = uz*ux*(1-np.cos(theta)) - uy*np.sin(theta)
    rot_mat[2][1] = uz*uy*(1-np.cos(theta)) + ux*np.sin(theta)
    rot_mat[2][2] = np.cos(theta) + uz**2*(1-np.cos(theta))
    
    return np.matmul(rot_mat, coord.T).T

def check_Cn(coord):
    '''
    Find the highest Cn (from C8 - C3) in the molecule along x, y or z axis, excluding C2 and C1.

    Args:
        coord: coordinates of the molecule

    Returns:
        order of rotation, n, if the molecule has Cn symmetry, None otherwise.
    '''

    for i in range(3, 9).__reversed__():
        rot_z = sort_coords(thre_cut(Cn_rot(coord, i, 'z'), 4))
        rot_y = sort_coords(thre_cut(Cn_rot(coord, i, 'y'), 4))
        rot_x = sort_coords(thre_cut(Cn_rot(coord, i, 'x'), 4))
        coord = sort_coords(thre_cut(coord, 4))
        # print(f"rot_z: {rot_z}")
        # print(f"rot_y: {rot_y}")
        # print(f"rot_x: {rot_x}")
        # print(f"coord: {coord}")

        if np.array_equal(rot_z, coord):
            print(f"C{i} rotation about z-axis found in the molecule.")
            return i
        elif np.array_equal(rot_y, coord):
            print(f"C{i} rotation about y-axis found in the molecule.")
            return i
        elif np.array_equal(rot_x, coord):
            print(f"C{i} rotation about x-axis found in the molecule.")
            return i
        else:
            continue
    
    print("No Cn (3-8) rotation axis along x, y, z found in the molecule.")
    return 0
    
def check_Cn_extd(coord):
    '''
    Find the highest Cn (from C8 - C3) in the molecule along some axis, excluding C2 and C1.

    the axis are:
    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]
    [ 1, -1, 0], [ 1, 0, -1], [0,  1, -1]
    [-1, -1, 0], [-1, 0, -1], [0, -1, -1]
    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 
    [-1, -1, 1], [-1, 1, -1], [1, -1, -1], 
    [-1, -1, -1]

        Args:
        coord: coordinates of the molecule

    Returns:
        order of rotation, n, if the molecule has Cn symmetry, None otherwise.
    '''

    for i in range(3, 9).__reversed__():
        for axis in extend_axis:
            axis = np.array(axis)/np.linalg.norm(axis)
            rot = sort_coords(thre_cut(Cn_rot_axis(coord, i, axis), 4))
            coord = sort_coords(thre_cut(coord, 4))
            if np.array_equal(rot, coord):
                print(f"C{i} rotation about {axis} found in the molecule.")
                return i
            else:
                continue
    print("No Cn (3-8) rotation axis found in the molecule.")
    return 0

# Improper Rotation

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
    irot_mat=np.matmul(rot_mat, ref_mat)
    return np.matmul(irot_mat, coord.T).T
    
def Sn_rot_axis(coord, n, axis):
    '''
    Rotate and reflects the coordinates by Sn along a given axis unit vector.

    Args:
        coord: coordinates of the molecule
        n: order of rotation, 2, 3, 4, 5, 6, ...
        axis: unit vector of axis of rotation, [ux, uy, uz], where ux^2 + uy^2 + uz^2 = 1

    Rotational matrix:
    [[cos(theta) + ux^2*(1-cos(theta)),    ux*uy*(1-cos(theta)) - uz*sin(theta), ux*uz*(1-cos(theta)) + uy*sin(theta)],
    [uy*ux*(1-cos(theta)) + uz*sin(theta), cos(theta) + uy^2*(1-cos(theta)),     uy*uz*(1-cos(theta)) - ux*sin(theta)],
    [uz*ux*(1-cos(theta)) - uy*sin(theta), uz*uy*(1-cos(theta)) + ux*sin(theta), cos(theta) + uz^2*(1-cos(theta))]]

    where theta = 2*pi/n

    Reflection matrix:
    [[1 - 2*ux^2, -2*ux*uy, -2*ux*uz],
    [-2*uy*ux, 1 - 2*uy^2, -2*uy*uz],
    [-2*uz*ux, -2*uz*uy, 1 - 2*uz^2]]

    '''
    rot_mat = np.zeros((3,3))
    ref_mat = np.zeros((3,3))
    theta = 2*np.pi/n
    ux, uy, uz = axis
    rot_mat[0][0] = np.cos(theta) + ux**2*(1-np.cos(theta))
    rot_mat[0][1] = ux*uy*(1-np.cos(theta)) - uz*np.sin(theta)
    rot_mat[0][2] = ux*uz*(1-np.cos(theta)) + uy*np.sin(theta)
    rot_mat[1][0] = uy*ux*(1-np.cos(theta)) + uz*np.sin(theta)
    rot_mat[1][1] = np.cos(theta) + uy**2*(1-np.cos(theta))
    rot_mat[1][2] = uy*uz*(1-np.cos(theta)) - ux*np.sin(theta)
    rot_mat[2][0] = uz*ux*(1-np.cos(theta)) - uy*np.sin(theta)
    rot_mat[2][1] = uz*uy*(1-np.cos(theta)) + ux*np.sin(theta)
    rot_mat[2][2] = np.cos(theta) + uz**2*(1-np.cos(theta))

    ref_mat[0][0] = 1 - 2*ux**2
    ref_mat[0][1] = -2*ux*uy
    ref_mat[0][2] = -2*ux*uz
    ref_mat[1][0] = -2*uy*ux
    ref_mat[1][1] = 1 - 2*uy**2
    ref_mat[1][2] = -2*uy*uz
    ref_mat[2][0] = -2*uz*ux
    ref_mat[2][1] = -2*uz*uy
    ref_mat[2][2] = 1 - 2*uz**2

    irot_mat=np.matmul(rot_mat, ref_mat)
    return np.matmul(irot_mat, coord.T).T
        

def check_Sn(coord):
    '''
    Find the highest Sn (from S10 - S3) in the molecule along x, y or z axis, excluding S2 and S1.

    Args:
        coord: coordinates of the molecule

    Returns:
        order of rotation, n, if the molecule has Sn symmetry, None otherwise.
    '''

    for i in range(3, 11).__reversed__():
        irot_z = sort_coords(thre_cut(Sn_rot(coord, i, 'z'), 4))
        irot_y = sort_coords(thre_cut(Sn_rot(coord, i, 'y'), 4))
        irot_x = sort_coords(thre_cut(Sn_rot(coord, i, 'x'), 4))
        coord = sort_coords(thre_cut(coord, 4))
        # print(f"rot_z: {rot_z}")
        # print(f"rot_y: {rot_y}")
        # print(f"rot_x: {rot_x}")
        # print(f"coord: {coord}")

        if np.array_equal(irot_z, coord):
            print(f"S{i} improper rotation around z-axis found in the molecule.")
            return i
        elif np.array_equal(irot_y, coord):
            print(f"S{i} improper rotation around y-axis found in the molecule.")
            return i
        elif np.array_equal(irot_x, coord):
            print(f"S{i} improper rotation around x-axis found in the molecule.")
            return i
        else:
            continue
    
    print("No Sn (3-10) rotation axis along x, y, z found in the molecule.")
    return 0


def check_Sn_extd(coord):
    '''
    Find the highest Sn (from S10 - S3) in the molecule along some axis, excluding S2 and S1.

    the axis are:
    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]
    [ 1, -1, 0], [ 1, 0, -1], [0,  1, -1]
    [-1, -1, 0], [-1, 0, -1], [0, -1, -1]
    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 
    [-1, -1, 1], [-1, 1, -1], [1, -1, -1], 
    [-1, -1, -1]

        Args:
        coord: coordinates of the molecule

    Returns:
        order of rotation, n, if the molecule has Cn symmetry, None otherwise.
    '''

    for i in range(3, 11).__reversed__():
        for axis in extend_axis:
            axis = np.array(axis)/np.linalg.norm(axis)
            rot = sort_coords(thre_cut(Sn_rot_axis(coord, i, axis), 4))
            coord = sort_coords(thre_cut(coord, 4))
            if np.array_equal(rot, coord):
                print(f"S{i} improper rotation about {axis} found in the molecule.")
                return i
            else:
                continue
    print("No Sn (3-10) rotation axis found in the molecule.")
    return 0


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
    
    return np.matmul(ref_mat, coord.T).T

def reflection_plane(coord, vec):
    '''
    Reflect the coordinates about a plane with normal vector vec[ux, uy, uz].

    Standard reflection matrix:
    [[1 - 2*ux^2, -2*ux*uy, -2*ux*uz],
    [-2*uy*ux, 1 - 2*uy^2, -2*uy*uz],
    [-2*uz*ux, -2*uz*uy, 1 - 2*uz^2]]

    Args:
        coord: coordinates of the molecule
        vec: normal vector of the plane

    Returns:
        reflected coordinates of the molecule
    '''
    ref_mat = np.zeros((3,3))
    vec = np.array(vec)/np.linalg.norm(vec)
    for i in range(3):
        for j in range(3):
            ref_mat[i][j] = (delta(i,j) - 2*vec[i]*vec[j])

    return np.matmul(ref_mat, coord.T).T


def check_reflection(coord):
    '''
    Check the symmetry of the molecule using reflection.

    Args:
        coord: coordinates of the molecule

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
    coord = sort_coords(thre_cut(coord, 4))
    ref_xy = sort_coords(thre_cut(reflection(coord, 'xy'), 4))
    ref_yz = sort_coords(thre_cut(reflection(coord, 'yz'), 4))
    ref_xz = sort_coords(thre_cut(reflection(coord, 'xz'), 4))

    if np.array_equal(ref_xy, coord):
        print("Reflection about xy-plane found in the molecule.")
        return True
    elif np.array_equal(ref_yz, coord):
        print("Reflection about yz-plane found in the molecule.")
        return True
    elif np.array_equal(ref_xz, coord):
        print("Reflection about xz-plane found in the molecule.")
        return True
    else:
        print("No reflection mirror about xyz plane found in the molecule.")
        return False

def check_reflection_extd(coord):
    '''
    Check the symmetry of the molecule using reflection about somes planes.

    The unit vectors of the normal vectors of the planes are:

    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]
    [ 1, -1, 0], [ 1, 0, -1], [0,  1, -1]
    [-1, -1, 0], [-1, 0, -1], [0, -1, -1]
    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 
    [-1, -1, 1], [-1, 1, -1], [1, -1, -1], 
    [-1, -1, -1]
    
    Args:
        coord: coordinates of the molecule

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
    coord = sort_coords(thre_cut(coord, 4))
    for vec in extend_axis:
        ref = sort_coords(thre_cut(reflection_plane(coord, vec), 4))
        if np.array_equal(ref, coord):
            print(f"Reflection about plane with normal vector {vec} found in the molecule.")
            return True
        else:
            continue
    print("No reflection mirror found in the molecule.")
    return False

# Inversion

def inversion(coord):
    '''
    Invert the molecule

    Inversion matrix
    [[-1,  0,  0],
     [ 0, -1,  0],
     [ 0,  0, -1]]

    Args:
        coord: coordinates of the molecule

    Returns:
        inverted coordinates of the molecule
    '''
    inv_mat = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    inv_coord = np.matmul(inv_mat, coord.T).T

    return inv_coord


def check_inversion(coord):
    '''
    Check if there is an inversion center.

    Args:
        coord: coordinates of the molecule

    Returns:
        True if there is an inversion center, False otherwise.
    '''
    coord = sort_coords(thre_cut(coord, 4))
    inv_coord = sort_coords(thre_cut(inversion(coord), 4))

    if np.array_equal(inv_coord, coord):
        print("Inversion center found in the molecule.")
        return True
    else:
        print("No inversion center found in the molecule.")
        return False


# def invert(coord,atom):
#     xi=-coord[0]
#     yi=-coord[1]
#     zi=-coord[2]
#     xyzi=[float("{:.4f}".format(xi)), float("{:.4f}".format(yi)), float("{:.4f}".format(zi)),atom]
# #     return xyzi

# def test_inversion(cd):
#     atm_flg=[]
#     atm=[]
#     inv_atm=[]
#     final=0
#     for i in range (0,(cd.natm)):
#         atm_flg.insert(i,False)
#         atm.append(cd.coordinates[i])
#         atm.append(cd.atm_name[i])
#         #even positions stand for coordinates
#         #odd positions stand for atom name
#         if i%2!=0:
#             inv_atm.append(invert(atm[i-1],atm[i]))
    
#     for j in range (0,len(inv_atm)):
#         #check if it's possible to put j and k in the same for
#         for k in range (0,len(inv_atm)):
#             if k%2==0:
#                 if inv_atm[j][0:3]==atm[k].tolist():
#                     #add check the type of atom
#                     atm_flg[j]=True
   
#     for l in range(0,len(atm_flg)):
#         if atm_flg[l]==False:
#             final=final+1
    
#     if final==0:
#         print("There is an inversion center")
#         inv=True
#     else:
#         print("No inversion center")
#         inv=False

#     return 

