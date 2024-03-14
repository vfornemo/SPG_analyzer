# symmetry operation functions

from data import TOLERANCE
import math
from spg import Molecule


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


# Rotation
#for H2O molecule

# Rotation matrix function around the 3 main axis
def rotate_matrix_x (x, y, z, angle, x_shift=0, y_shift=0, z_shift=0, units="DEGREES"):
    # """
    # Rotates a point in the xy-plane counterclockwise through an angle about the origin
    # https://en.wikipedia.org/wiki/Rotation_matrix
    # :param x: x coordinate
    # :param y: y coordinate
    # :param z: z coordinate
    # :param x_shift: x-axis shift from origin (0, 0)
    # :param y_shift: y-axis shift from origin (0, 0)
    # :param z_shift: z-axis shift from origin (0, 0)
    # :param angle: The rotation angle in degrees
    # :param units: DEGREES (default) or RADIANS
    # :return: Tuple of rotated x, y, and z
    # """

    # Shift to origin (0,0)
    x = x - x_shift
    y = y - y_shift
    z = z - z_shift

    # Convert degrees to radians
    if units == "DEGREES":
        angle = math.radians(angle)

    # Rotation matrix multiplication to get rotated x, y & z
    xr = x + x_shift
    yr = (y * math.cos(angle)) - (z * math.sin(angle)) + y_shift
    zr = (y * math.sin(angle)) + (z * math.cos(angle)) + z_shift

    return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr))
def rotate_matrix_y (x, y, z, angle, x_shift=0, y_shift=0, z_shift=0, units="DEGREES"):
    # """
    # Rotates a point in the xy-plane counterclockwise through an angle about the origin
    # https://en.wikipedia.org/wiki/Rotation_matrix
    # :param x: x coordinate
    # :param y: y coordinate
    # :param z: z coordinate
    # :param x_shift: x-axis shift from origin (0, 0)
    # :param y_shift: y-axis shift from origin (0, 0)
    # :param z_shift: z-axis shift from origin (0, 0)
    # :param angle: The rotation angle in degrees
    # :param units: DEGREES (default) or RADIANS
    # :return: Tuple of rotated x, y, and z
    # """

    # Shift to origin (0,0)
    x = x - x_shift
    y = y - y_shift
    z = z - z_shift

    # Convert degrees to radians
    if units == "DEGREES":
        angle = math.radians(angle)

    # Rotation matrix multiplication to get rotated x, y & z
    xr = (x * math.cos(angle)) + (z * math.sin(angle)) + x_shift
    yr = y + y_shift
    zr = (-x * math.sin(angle)) + (z * math.cos(angle)) + z_shift

    return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr))
def rotate_matrix_z (x, y, z, angle, x_shift=0, y_shift=0, z_shift=0, units="DEGREES"):
    # """
    # Rotates a point in the xy-plane counterclockwise through an angle about the origin
    # https://en.wikipedia.org/wiki/Rotation_matrix
    # :param x: x coordinate
    # :param y: y coordinate
    # :param z: z coordinate
    # :param x_shift: x-axis shift from origin (0, 0)
    # :param y_shift: y-axis shift from origin (0, 0)
    # :param z_shift: z-axis shift from origin (0, 0)
    # :param angle: The rotation angle in degrees
    # :param units: DEGREES (default) or RADIANS
    # :return: Tuple of rotated x, y, and z
    # """

    # Shift to origin (0,0)
    x = x - x_shift
    y = y - y_shift
    z = z - z_shift

    # Convert degrees to radians
    if units == "DEGREES":
        angle = math.radians(angle)

    # Rotation matrix multiplication to get rotated x, y & z
    xr = (x * math.cos(angle)) - (y * math.sin(angle)) + x_shift
    yr = (x * math.sin(angle)) + (y * math.cos(angle)) + y_shift
    zr = z + z_shift

    return float("{:.4f}".format(xr)), float("{:.4f}".format(yr)), float("{:.4f}".format(zr))
#Testing C2 for water molecule

#1 - H1
#2 - O
#3 - H2

atom_type1='H'
atom_type2='O'
atom_type3='H'

x1=0.8110
y1=-0.4677
z1=-0.0000
x2=-0.0000
y2=0.0589
z2=0.0000
x3=-0.8110
y3=-0.4677
z3=-0.0000

#relate atom type to coordinates

def test_Cn(x1,y1,z1,x2,y2,z2,x3,y3,z3):

    angles=[180, 120, 90]
    higher=1
    for angle in angles:
        atom1=False
        atom2=False
        atom3=False
        #performing the rotation
        #next step - get the coordinates directly from the .mol file
        coord1_r=rotate_matrix_x(x1,y1,z1,angle)
        coord2_r=rotate_matrix_x(x2,y2,z2,angle)
        coord3_r=rotate_matrix_x(x3,y3,z3,angle)
        #rename the elements
        x1_r=coord1_r[0]
        y1_r=coord1_r[1]
        z1_r=coord1_r[2]

        x2_r=coord2_r[0]
        y2_r=coord2_r[1]
        z2_r=coord2_r[2]

        x3_r=coord3_r[0]
        y3_r=coord3_r[1]
        z3_r=coord3_r[2]
    #then, operate all coordinates at once
    #loop for
        #testing if this operation of symmetry exists
        #check if O is at the same place & if we find H where it used to be before
        if x2_r == x2 and y2_r == y2 and z2_r == z2:
            atom1=True
            #print("O is ok")
        # else:
        #     print("O is NOT ok")
        if x1_r == x3 and y1_r == y3 and z1_r == z3:
            atom2=True
            #print("H1 is ok")
        # else:
        #     print("H1 is NOT ok")
        if x3_r == x1 and y3_r == y1 and z3_r == z1:
            atom3=True
            #print("H2 is ok")
        # else:
        #     print("H2 is NOT ok")
        if atom1==True and atom2==True and atom3==True:
            print(f"There is a C{int(360/angle)} symmetry operation")
            main_axis='x'
            if int(360/angle)>higher:
                higher=int(360/angle)
    for angle in angles:
        atom1=False
        atom2=False
        atom3=False
        #performing the rotation
        #next step - get the coordinates directly from the .mol file
        coord1_r=rotate_matrix_y(x1,y1,z1,angle)
        coord2_r=rotate_matrix_y(x2,y2,z2,angle)
        coord3_r=rotate_matrix_y(x3,y3,z3,angle)
        #rename the elements
        x1_r=coord1_r[0]
        y1_r=coord1_r[1]
        z1_r=coord1_r[2]

        x2_r=coord2_r[0]
        y2_r=coord2_r[1]
        z2_r=coord2_r[2]

        x3_r=coord3_r[0]
        y3_r=coord3_r[1]
        z3_r=coord3_r[2]
    #then, operate all coordinates at once
    #loop for
        #testing if this operation of symmetry exists
        #check if O is at the same place & if we find H where it used to be before
        if x2_r == x2 and y2_r == y2 and z2_r == z2:
            atom1=True
            #print("O is ok")
        # else:
        #     print("O is NOT ok")
        if x1_r == x3 and y1_r == y3 and z1_r == z3:
            atom2=True
            #print("H1 is ok")
        # else:
        #     print("H1 is NOT ok")
        if x3_r == x1 and y3_r == y1 and z3_r == z1:
            atom3=True
            #print("H2 is ok")
        # else:
        #     print("H2 is NOT ok")
        if atom1==True and atom2==True and atom3==True:
            print(f"There is a C{int(360/angle)} symmetry operation")
            main_axis='y'
            if int(360/angle)>higher:
                higher=int(360/angle)
    for angle in angles:
        atom1=False
        atom2=False
        atom3=False
        #performing the rotation
        #next step - get the coordinates directly from the .mol file
        coord1_r=rotate_matrix_z(x1,y1,z1,angle)
        coord2_r=rotate_matrix_z(x2,y2,z2,angle)
        coord3_r=rotate_matrix_z(x3,y3,z3,angle)
        #rename the elements
        x1_r=coord1_r[0]
        y1_r=coord1_r[1]
        z1_r=coord1_r[2]

        x2_r=coord2_r[0]
        y2_r=coord2_r[1]
        z2_r=coord2_r[2]

        x3_r=coord3_r[0]
        y3_r=coord3_r[1]
        z3_r=coord3_r[2]
    #then, operate all coordinates at once
    #loop for
        #testing if this operation of symmetry exists
        #check if O is at the same place & if we find H where it used to be before
        if x2_r == x2 and y2_r == y2 and z2_r == z2:
            atom1=True
            #print("O is ok")
        # else:
        #     print("O is NOT ok")
        if x1_r == x3 and y1_r == y3 and z1_r == z3:
            atom2=True
            #print("H1 is ok")
        # else:
        #     print("H1 is NOT ok")
        if x3_r == x1 and y3_r == y1 and z3_r == z1:
            atom3=True
            #print("H2 is ok")
        # else:
        #     print("H2 is NOT ok")
        if atom1==True and atom2==True and atom3==True:
            print(f"There is a C{int(360/angle)} symmetry operation")
            main_axis='z'
            if int(360/angle)>higher:
                higher=int(360/angle)
    if higher!=1:
        print(f"The higher rotation symmetry operation is C{higher} around {main_axis}")
    else:
        print("There is no rotation symmetry operation")
    
    return

test_Cn(x1,y1,z1,x2,y2,z2,x3,y3,z3)

#to do:
#get the coordinates directly from the .mol file
#operate all coordinates at once (no 1, 2, 3... manually)
#attribute atom types to the coordinates (no 1, 2, 3...), so it's possible to check if there is the same atom at the desired position - generalization

