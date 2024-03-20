from calendar import c
from math import e
import numpy as np
# symmetry operation functions

from utils import delta, sort_coords, sort_coords_idx, thre_cut, compare_coords
from data import extend_axis, full_axis

# reference: http://www.pci.tu-bs.de/aggericke/PC4e/Kap_IV/Matrix_Symm_Op.htm
# Proper Rotation

# Cn

# reference: https://en.wikipedia.org/wiki/Rotation_matrix
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

# reference: https://en.wikipedia.org/wiki/Rotation_matrix
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

def check_Cn(atm_list, coord, return_axis=False, mode='medium', enhanced=False):
    '''
    Find the highest Cn (from C8 - C2) in the molecule along x, y or z axis, excluding C1.

    Args:
        atm_list: list of atoms in the molecule
        coord: coordinates of the molecule
        return_axis: if True, return the axis of rotation
        mode: tolerance mode, from 'ultra_loose' to 'very_tight', default is 'medium'
        enhanced: if True, the comparison will enable coords intersection function

    Returns:
        order of rotation, n, if the molecule has Cn symmetry, 0 otherwise.
        vector of the Cn axis, [ux, uy, uz]
    '''

    for i in range(2, 9).__reversed__():
        rot_z = Cn_rot(coord, i, 'z')
        rot_y = Cn_rot(coord, i, 'y')
        rot_x = Cn_rot(coord, i, 'x')

        # print(f"Check C{i} rotation about x, y, z axis.")
        # rot_z1 = thre_cut(rot_z, 2)
        # rot_y1 = thre_cut(rot_y, 2)
        # rot_x1 = thre_cut(rot_x, 2)
        # coord1 = thre_cut(coord, 2)
        # print(f"rot_z:", sort_coords(rot_z1, sort_coords_idx(rot_z1)))
        # print(f"rot_y:", sort_coords(rot_y1, sort_coords_idx(rot_y1)))
        # print(f"rot_x:", sort_coords(rot_x1, sort_coords_idx(rot_x1)))
        # print(f"coord:", sort_coords(coord1, sort_coords_idx(coord1)))

        if compare_coords(atm_list, rot_z, coord, mode, enhanced):
            print(f"C{i} rotation about z-axis found in the molecule.")
            if return_axis:
                return i, [0, 0, 1]
            else:
                return i
        elif compare_coords(atm_list, rot_y, coord, mode, enhanced):
            print(f"C{i} rotation about y-axis found in the molecule.")
            if return_axis:
                return i, [0, 1, 0]
            else:
                return i
        elif compare_coords(atm_list, rot_x, coord, mode, enhanced):
            print(f"C{i} rotation about x-axis found in the molecule.")
            if return_axis:
                return i, [1, 0, 0]
            else:
                return i
        else:
            continue
    
    print("No Cn (2-8) rotation axis along x, y, z found in the molecule.")

    if return_axis:
        return 0, [0, 0, 0]
    else:
        return 0
    
def check_Cn_extd(atm_list, coord, return_axis=False, mode='medium', enhanced=False):
    '''
    Find the highest Cn (from C8 - C2) in the molecule along some axis, and C1.

    the axis are:
    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]

    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 


    Args:
    coord: coordinates of the molecule
    return_axis: if True, return the axis of rotation

    Returns:
        order of rotation, n, if the molecule has Cn symmetry, 0 otherwise.
        vector of the Cn axis, [ux, uy, uz]
    '''

    for i in range(2, 9).__reversed__():
        for axis in extend_axis:
            axis = np.array(axis)/np.linalg.norm(axis)
            rot = Cn_rot_axis(coord, i, axis)

            if compare_coords(atm_list, rot, coord, mode, enhanced):
                print(f"C{i} rotation about {axis} found in the molecule.")
                if return_axis:
                    return i, axis
                else:
                    return i
            else:
                continue

    print("No Cn (2-8) rotation axis found in the molecule.")
    if return_axis:
        return 0, [0, 0, 0]
    else:
        return 0

def check_C2_perp_Cn(atm_list, coord, Cn, mode='medium', enhanced=False):
    '''
    Check how many C2 axis that are perpendicular to the Cn axis in the molecule.

    Args:
        coord: coordinates of the molecule
        Cn: vector of Cn axis

    Returns:
        number of C2 axis, n, if the molecule has C2 symmetry, 0 otherwise.
    '''

    if Cn == [0, 0, 0]:
        print("No Cn axis in the molecule.")
        return 0
    
    nC2 = 0

    for axis in full_axis:
        if np.dot(Cn, axis) == 0:
            axis = np.array(axis)/np.linalg.norm(axis)
            rot = Cn_rot_axis(coord, 2, axis)
            if compare_coords(atm_list, rot, coord, mode, enhanced):
                nC2 += 1

    print(f"{nC2} C2 axis perpendicular to {Cn} found in the molecule.")
    return nC2

def check_C2_full(atm_list, coord, mode='medium'):
    '''
    Check how many C2 axis that are perpendicular to the Cn axis in the molecule.

    Args:
        coord: coordinates of the molecule
        Cn: vector of Cn axis

    Returns:
        number of C2 axis, n, if the molecule has C2 symmetry, 0 otherwise.
    '''
    nC2 = 0

    for axis in full_axis:
        axis = np.array(axis)/np.linalg.norm(axis)
        rot = Cn_rot_axis(coord, 2, axis)
        if compare_coords(atm_list, rot, coord, mode):
            nC2 += 1

    print(f"{nC2} C2 axis found in the molecule.")
    return nC2

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
        

def check_Sn(atm_list, coord, check_from=8, mode='medium'):
    '''
    Find the highest Sn (from S{check_from} - S3) in the molecule along x, y or z axis, excluding S2 and S1.

    Args:
        coord: coordinates of the molecule
        check_from: the highest order of rotation to check

    Returns:
        order of rotation, n, if the molecule has Sn symmetry, 0 otherwise.
    '''

    for i in range(3, check_from+1).__reversed__():
        irot_z = Sn_rot(coord, i, 'z')
        irot_y = Sn_rot(coord, i, 'y')
        irot_x = Sn_rot(coord, i, 'x')
        # print(f"rot_z: {rot_z}")
        # print(f"rot_y: {rot_y}")
        # print(f"rot_x: {rot_x}")
        # print(f"coord: {coord}")

        if compare_coords(atm_list, irot_z, coord, mode):
            print(f"S{i} improper rotation around z-axis found in the molecule.")
            return i
        elif compare_coords(atm_list, irot_y, coord, mode):
            print(f"S{i} improper rotation around y-axis found in the molecule.")
            return i
        elif compare_coords(atm_list, irot_x, coord, mode):
            print(f"S{i} improper rotation around x-axis found in the molecule.")
            return i
        else:
            continue
    
    print(f"No Sn (3-{check_from}) rotation axis along x, y, z found in the molecule.")
    return 0


def check_Sn_extd(atm_list, coord, check_from=8, mode='medium'):
    '''
    Find the highest Sn (from S{check_from} - S3) in the molecule along some axis, excluding S2 and S1.

    the axis are:
    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]

    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 


    Args:
        coord: coordinates of the molecule
        check_from: the highest order of rotation to check

    Returns:
        order of rotation, n, if the molecule has Sn symmetry, 0 otherwise.
    '''

    for i in range(3, check_from+1).__reversed__():
        for axis in extend_axis:
            axis = np.array(axis)/np.linalg.norm(axis)
            rot = Sn_rot_axis(coord, i, axis)
            if compare_coords(atm_list, rot, coord, mode):
                print(f"S{i} improper rotation about {axis} found in the molecule.")
                return i
            else:
                continue
    print(f"No Sn (3-{check_from}) rotation axis found in the molecule.")
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

# reference: https://en.wikipedia.org/wiki/Transformation_matrix
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


def check_reflection(atm_list, coord, mode='medium'):
    '''
    Check the symmetry of the molecule using reflection.

    Args:
        coord: coordinates of the molecule

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
    ref_xy = reflection(coord, 'xy')
    ref_yz = reflection(coord, 'yz')
    ref_xz = reflection(coord, 'xz')

    if compare_coords(atm_list, ref_xy, coord, mode):
        print("Reflection about xy-plane found in the molecule.")
        return True
    elif compare_coords(atm_list, ref_yz, coord, mode):
        print("Reflection about yz-plane found in the molecule.")
        return True
    elif compare_coords(atm_list, ref_xz, coord, mode):
        print("Reflection about xz-plane found in the molecule.")
        return True
    else:
        print("No reflection mirror about xyz plane found in the molecule.")
        return False

def check_reflection_extd(atm_list, coord, mode='medium'):
    '''
    Check the symmetry of the molecule using reflection about somes planes.

    The unit vectors of the normal vectors of the planes are:

    1. axis between x, y, z
    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
    [-1,  1, 0], [-1, 0,  1], [0, -1,  1]

    2. axis among xyz axis
    [1, 1, 1], 
    [-1, 1, 1], [1, -1, 1], [1, 1, -1], 

    
    Args:
        coord: coordinates of the molecule

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
    for vec in extend_axis:
        ref = reflection_plane(coord, vec)
        if compare_coords(atm_list, ref, coord, mode):
            print(f"Reflection about plane with normal vector {vec} found in the molecule.")
            return True
        else:
            continue
    print("No reflection mirror found in the molecule.")
    return False

def check_reflection_h(atm_list, coord, Cn, mode='medium'):
    '''
    Check the symmetry of the molecule using reflection about a plane perpendicular to the Cn axis (sigma_h).

    Args:
        coord: coordinates of the molecule
        Cn: vector of Cn axis

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
        
    if Cn == [0, 0, 0]:
        print("No Cn axis in the molecule.")
        return  

    ref = reflection_plane(coord, Cn)
    if compare_coords(atm_list, ref, coord, mode):
        print(f"Reflection about plane perpendicular to {Cn} found in the molecule.")
        return True
    else:
        print("No reflection mirror h found in the molecule.")
        return False	  

def check_reflection_v(atm_list, coord, Cn, mode='medium'):
    '''
    Check the symmetry of the molecule using reflection about a plane including the Cn axis (sigma_v).
    Three planes, xy, yz, xz, are considered. The ones that are not perpendicular to the Cn axis are considered.

    Args:
        coord: coordinates of the molecule
        Cn: vector of Cn axis

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''
        
    if Cn == [0, 0, 0]:
        print("No Cn axis in the molecule.")
        return  
    
    plane_vec = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    for vec in plane_vec:
        if np.dot(Cn, vec) == 0:
            ref = reflection_plane(coord, vec)
            if compare_coords(atm_list, ref, coord, mode):
                print(f"Reflection about plane including {vec} found in the molecule.")
                return True
            else:
                continue
        else:
            continue
        
    print("No reflection mirror v found in the molecule.")
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


def check_inversion(atm_list, coord, mode='medium'):
    '''
    Check if there is an inversion center.

    Args:
        coord: coordinates of the molecule

    Returns:
        True if there is an inversion center, False otherwise.
    '''
    inv_coord = inversion(coord)

    if compare_coords(atm_list, inv_coord, coord, mode):
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

