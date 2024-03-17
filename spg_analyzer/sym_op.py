from calendar import c
import numpy as np
# symmetry operation functions

from utils import sort_coords, thre_cut


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

def check_Cn(coord):
    '''
    Check the symmetry of the molecule using Cn (from C8 - C3), excluding C2 and C1.

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
            print(f"C{i} rotation about z-axis is a symmetry operation for the molecule.")
            return i
        elif np.array_equal(rot_y, coord):
            print(f"C{i} rotation about y-axis is a symmetry operation for the molecule.")
            return i
        elif np.array_equal(rot_x, coord):
            print(f"C{i} rotation about x-axis is a symmetry operation for the molecule.")
            return i
        else:
            continue
    
    print("No Cn (3-8) rotation is a symmetry operation for the molecule.")
    return None
    

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
        elif np.array_equal(irot_y, coord):
            print(f"S{i} improper rotation around y-axis is a symmetry operation for the molecule.")
            return i
        elif np.array_equal(irot_x, coord):
            print(f"S{i} improper rotation around x-axis is a symmetry operation for the molecule.")
            return i
        else:
            continue
    
    print("No Sn (3-8) rotation is a symmetry operation for the molecule.")
    return None

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


def check_reflection(coord):
    '''
    Check the symmetry of the molecule using reflection.

    Args:
        coord: coordinates of the molecule

    Returns:
        True if the molecule has reflection symmetry, False otherwise.
    '''

    ref_xy = thre_cut(reflection(coord, 'xy'))
    ref_yz = thre_cut(reflection(coord, 'yz'))
    ref_xz = thre_cut(reflection(coord, 'xz'))

    if np.array_equal(sort_coords(ref_xy), sort_coords(coord)):
        print("Reflection about xy-plane is a symmetry operation for the molecule.")
        return True
    elif np.array_equal(sort_coords(ref_yz), sort_coords(coord)):
        print("Reflection about yz-plane is a symmetry operation for the molecule.")
        return True
    elif np.array_equal(sort_coords(ref_xz), sort_coords(coord)):
        print("Reflection about xz-plane is a symmetry operation for the molecule.")
        return True
    else:
        print("Reflection is not a symmetry operation for the molecule.")
        return False

# Inversion

def invert(coord,atom):
    xi=-coord[0]
    yi=-coord[1]
    zi=-coord[2]
    xyzi=[float("{:.4f}".format(xi)), float("{:.4f}".format(yi)), float("{:.4f}".format(zi)),atom]
    return xyzi

def test_inversion(cd):
    atm_flg=[]
    atm=[]
    inv_atm=[]
    final=0
    for i in range (0,(cd.natm)):
        atm_flg.insert(i,False)
        atm.append(cd.coordinates[i])
        atm.append(cd.atm_name[i])
        #even positions stand for coordinates
        #odd positions stand for atom name
        if i%2!=0:
            inv_atm.append(invert(atm[i-1],atm[i]))
    
    for j in range (0,len(inv_atm)):
        #check if it's possible to put j and k in the same for
        for k in range (0,len(inv_atm)):
            if k%2==0:
                if inv_atm[j][0:3]==atm[k].tolist():
                    #add check the type of atom
                    atm_flg[j]=True
   
    for l in range(0,len(atm_flg)):
        if atm_flg[l]==False:
            final=final+1
    
    if final==0:
        print("There is an inversion center")
        inv=True
    else:
        print("No inversion center")
        inv=False

    return 

