from calendar import c
import numpy as np
# symmetry operation functions

from utils import sort_coords, thre_cut


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


mol1 = Molecule([['H', [0.0, 0.0, 0.0]], 
   ['H', [0.0, 0.0, 1.0]],
   ['O', [0.0, 1.0, 0.0]]])

mol1.build()

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

test_inversion(mol1)


# # Rotation


# # Rotation matrix function around the 3 main axis
# def rotate_matrix_x (x, y, z, angle, atom):

#     # Shift to origin (0,0)
#     x = x - x_shift
#     y = y - y_shift
#     z = z - z_shift

#     # Rotation matrix multiplication to get rotated x, y & z
#     xr = x + x_shift
#     yr = (y * math.cos(angle)) - (z * math.sin(angle)) + y_shift
#     zr = (y * math.sin(angle)) + (z * math.cos(angle)) + z_shift

#     return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr)), atom
# def rotate_matrix_y (x, y, z, angle, atom, units="DEGREES"):

#     # Shift to origin (0,0)
#     x = x - x_shift
#     y = y - y_shift
#     z = z - z_shift

#     # Convert degrees to radians
#     if units == "DEGREES":
#         angle = math.radians(angle)

#     # Rotation matrix multiplication to get rotated x, y & z
#     xr = (x * math.cos(angle)) + (z * math.sin(angle)) + x_shift
#     yr = y + y_shift
#     zr = (-x * math.sin(angle)) + (z * math.cos(angle)) + z_shift

#     return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr)), atom
# def rotate_matrix_z (x, y, z, angle, atom, units="DEGREES"):

#     # Shift to origin (0,0)
#     x = x - x_shift
#     y = y - y_shift
#     z = z - z_shift

#     # Convert degrees to radians
#     if units == "DEGREES":
#         angle = math.radians(angle)

#     # Rotation matrix multiplication to get rotated x, y & z
#     xr = (x * math.cos(angle)) - (y * math.sin(angle)) + x_shift
#     yr = (x * math.sin(angle)) + (y * math.cos(angle)) + y_shift
#     zr = z + z_shift

#     return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr)), atom

# mol1 = Molecule([['H', [0.0, 0.0, 0.0]], 
#    ['H', [0.0, 0.0, 1.0]],
#    ['O', [0.0, 1.0, 0.0]]])

# mol1.build()

# def test_Cn(cd):

#     angles=[2*math.pi/1,  360/2, 360/3, 360/4, 360/5, 360/6, 360/7, 360/8]
#     higher=1
#     atm_flgx=[]
#     atm_flgy=[]
#     atm_flgz=[]
#     atm=[]
#     rtt_atmx=[]
#     rtt_atmy=[]
#     rtt_atmz=[]
#     finalx=0
#     finaly=0
#     finalz=0

#     for i in angles:
#         for at in range (0,(cd.natm)):
#             atm_flg.insert(at,False)
#             atm.append(cd.coordinates[at])
#             atm.append(cd.atm_name[at])
#             #even positions stand for coordinates
#             #odd positions stand for atom name
#             if at%2!=0:
#                 rtt_atmx.append(rotate_matrix_x(atm[at-1],angles[at], atm[at]))
#                 rtt_atmy.append(rotate_matrix_y(atm[at-1],angles[at], atm[at]))
#                 rtt_atmz.append(rotate_matrix_z(atm[at-1],angles[at], atm[at]))
    
#     for j in range (0,len(rtt_atmx)):
#         #check if it's possible to put j and k in the same for
#         for k in range (0,len(rtt_atmx)):
#             if k%2==0:
#                 if rtt_atmx[j][0:3]==atm[k].tolist() and rtt_atmx[j][3]==atm[k+1]:
#                     atm_flgx[j]=True
#                 if rtt_atmy[j][0:3]==atm[k].tolist() and rtt_atmy[j][3]==atm[k+1]:
#                     atm_flgy[j]=True
#                 if rtt_atmz[j][0:3]==atm[k].tolist() and rtt_atmz[j][3]==atm[k+1]:
#                     atm_flgz[j]=True
    
#     for l in range(0,len(atm_flgx)):
#         if atm_flgx[l]==False:
#             finalx=finalx+1
#         if atm_flgy[l]==False:
#             finaly=finaly+1
#         if atm_flgz[l]==False:
#             finalz=finalz+1
    
#     for m in angles:
#         if finalx==0:
#             print(f"There is a C{int((m+1)*angles[m])} symmetry operation")
#             #main_axis.append('x')
#             rot_x.insert(m,int((m+1)*angles[m]))
#         if finaly==0:
#             print(f"There is a C{int((m+1)*angles[m])} symmetry operation")
#             #main_axis.append('y')
#             rot_y.insert(m,int((m+1)*angles[m]))
#         if finalz==0:
#             print(f"There is a C{int((m+1)*angles[m])} symmetry operation")
#             #main_axis.append('z')
#             rot_z.insert(m,int((m+1)*angles[m]))
#             if int((m+1)*angles[m])>higher:
#                 higher=int((m+1)*angles[m])

#     max_x=max(rot_x)
#     max_y=max(rot_y)
#     max_z=max(rot_z)
#     if max_x>max_y and max_x>max_z:
#         main_axis='x'
#         higher=max_x
#     elif max_y>max_x and max_y>max_z:
#         main_axis='y'
#         higher=max_y
#     elif max_z>max_x and max_z>max_y:
#         main_axis='z'
#         higher=max_z
    
#     if higher!=1:
#         print(f"The higher rotation symmetry operation is C{higher} around {main_axis}")
#     else:
#         print("There is no rotation symmetry operation")
    
#     return

# test_Cn(mol1)


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