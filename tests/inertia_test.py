import numpy as np


""" 
Dictionary of Periodic Table.
The keys are the symbols of elements, 
The values are lists of atomic numbers and masses (**INCOMPLETE**!).
""" 
dic_periodic_table = {
"H" :[1  , 1.00797],
"He":[2  , 4.0026],
"Li":[3  , 6.939],
"Be":[4  ,],
"B" :[5  ,],
"C" :[6  , 12.01115],
"N" :[7  , 14.0067],
"O" :[8  , 15.9994],
"F" :[9  ,],
"Ne":[10 ,],
"Na":[11 ,],
"Mg":[12 ,],
"Al":[13 ,],
"Si":[14 ,],
"P" :[15 ,],
"S" :[16 ,],
"Cl":[17 ,  34.9689],
"Ar":[18 ,],
"K" :[19 ,],
"Ca":[20 ,],
"Sc":[21 ,],
"Ti":[22 ,],
"V" :[23 ,],
"Cr":[24 ,],
"Mn":[25 ,],
"Fe":[26 ,],
"Co":[27 ,],
"Ni":[28 ,],
"Cu":[29 ,],
"Zn":[30 ,],
"Ga":[31 ,],
"Ge":[32 ,],
"As":[33 ,],
"Se":[34 ,],
"Br":[35 ,],
"Kr":[36 ,],
"Rb":[37 ,],
"Sr":[38 ,],
"Y" :[39 ,],
"Zr":[40 ,],
"Nb":[41 ,],
"Mo":[42 ,],
"Tc":[43 ,],
"Ru":[44 ,],
"Rh":[45 ,],
"Pd":[46 ,],
"Ag":[47 ,],
"Cd":[48 ,],
"In":[49 ,],
"Sn":[50 ,],
"Sb":[51 ,],
"Te":[52 ,],
"I" :[53 ,],
"Xe":[54 ,],
"Cs":[55 ,],
"Ba":[56 ,],
"La":[57 ,],
"Ce":[58 ,],
"Pr":[59 ,],
"Nd":[60 ,],
"Pm":[61 ,],
"Sm":[62 ,],
"Eu":[63 ,],
"Gd":[64 ,],
"Tb":[65 ,],
"Dy":[66 ,],
"Ho":[67 ,],
"Er":[68 ,],
"Tm":[69 ,],
"Yb":[70 ,],
"Lu":[71 ,],
"Hf":[72 ,],
"Ta":[73 ,],
"W" :[74 ,],
"Re":[75 ,],
"Os":[76 ,],
"Ir":[77 ,],
"Pt":[78 ,],
"Au":[79 ,],
"Hg":[80 ,],
"Tl":[81 ,],
"Pb":[82 ,],
"Bi":[83 ,],
"Po":[84 ,],
"At":[85 ,],
"Rn":[86 ,],
"Fr":[87 ,],
"Ra":[88 ,],
"Ac":[89 ,],
"Th":[90 ,],
"Pa":[91 ,],
"U" :[92 ,],
"Np":[93 ,],
"Pu":[94 ,],
"Am":[95 ,],
"Cm":[96 ,],
"Bk":[97 ,],
"Cf":[98 ,],
"Es":[99 ,],
"Fm":[100,],
"Md":[101,],
"No":[102,],
"Lr":[103,],
"Rf":[104,],
"Db":[105,],
"Sg":[106,],
"Bh":[107,],
"Hs":[108,],
"Mt":[109,],
"Ds":[110,],
"Rg":[111,],
"Cn":[112,],
"Nh":[113,],
"Fl":[114,],
"Mc":[115,],
"Lv":[116,],
"Ts":[117,],
"Og":[118,]}


def symbol_to_mass(coord):
    """ This function replace the element symbols i.e. the first
    column of Corrdinates to the corresponding mass. 

    It returns the new Corrdinates.
    
    ------------------------
    For instance, the input is
    H2O = [["O",  0.0,  0.0,  0.0],
           ["H",  0.7, -0.7,  0.0],
           ["H", -0.7, -0.7,  0.0]]
    ==>>
    output = [[15.9994    0.    0.0  0.  ]
              [1.00797    0.7  -0.7  0.  ]
              [1.00797   -0.7  -0.7  0.  ]]
    """
    for i in coord:
        if i[0] in dic_periodic_table.keys():
            i[0] = dic_periodic_table[i[0]][1]
    return np.asfarray(coord)
    

def calc_center_of_mass(data_array):
    """  This function calculates the center of mass of a given data array.
    The format of the data array is a two dimensional array, which contains
    the x, y, z Cartesion coordinates and the corresponding mass values for
    each x, y, z coordinates.

    It returns a 1-D data array.
    """
    return np.average(data_array[:,1:], axis=0, weights=data_array[:,0])


def calc_inertia_tensor(new_coord, mass):
    """  This function calculates the Elements of inertia tensor for the
    center-moved coordinates.
    
    The structure of the array is a two dimensional array, which contains
    the mass of elements (the first column) and their corresponded x, y, z
    Cartesion coordinates.

    It returns an inertia tensor.
    """
    I_xx = (mass * np.sum(np.square(new_coord[:,1:3:1]),axis=1)).sum()
    I_yy = (mass * np.sum(np.square(new_coord[:,0:3:2]),axis=1)).sum()
    I_zz = (mass * np.sum(np.square(new_coord[:,0:2:1]),axis=1)).sum()
    I_xy = (-1 * mass * np.prod(new_coord[:,0:2:1],axis=1)).sum()
    I_yz = (-1 * mass * np.prod(new_coord[:,1:3:1],axis=1)).sum()
    I_xz = (-1 * mass * np.prod(new_coord[:,0:3:2],axis=1)).sum()
    I = np.array([[I_xx, I_xy, I_xz],
		  [I_xy, I_yy, I_yz],
		  [I_xz, I_yz, I_zz]])
    return I


def find_principal_axes(data_array):
    """  This function finds the principal axes (I_a, I_b, I_c) by using
    diagonalizing the inertia tensor.

    It returns a diagonalized inertia tensor.
    """
    eig_val, eig_vec = np.linalg.eigh(data_array)
    print("Eigenvalues = \n{:}\n".format(eig_val))
    print("Eigenvectors = \n{:}\n".format(eig_vec))

    diag = np.dot(np.dot(np.linalg.inv(eig_vec), data_array), eig_vec)
    return np.around(diag, decimals=4)


def main(data_array):
    """ 
    workflow
    """
    print("Input coordinates = \n{:}\n".format(data_array))

    data_array = symbol_to_mass(data_array)
    CoM_coord = calc_center_of_mass(data_array)
    print("Center of Mass = \n{:}\n".format(CoM_coord))
    
    new_coord = data_array[:,1:] - CoM_coord
    print("shifted coordinates = \n{:}\n".format(np.column_stack((data_array[:,0],
                                                                  new_coord))))
    I = calc_inertia_tensor(new_coord, data_array[:,0])
    D = find_principal_axes(I)
    return I, D


###############################################################################
############################ some examples for test  ##########################
#### NH3
##a = [ 
##["N",  0     ,   0     ,   0],
##["H",  0     ,  -0.9377,  -0.3816],
##["H",  0.8121,   0.4689,  -0.3816],
##["H", -0.8121,   0.4689,  -0.3816]]

### CO2
##a = [
##["C",     -0.39214687,   -0.78171310,    0.01418586],
##["O",     -1.65054687,   -0.78171310,    0.01418586],
##["O",      0.86625313,   -0.78171310,    0.01418586]]

#### CH4
##a = [
##["C",        0,        0,        0],
##["H",   0.6276,   0.6276,   0.6276],
##["H",   0.6276,  -0.6276,  -0.6276],
##["H",  -0.6276,   0.6276,  -0.6276],
##["H",  -0.6276,  -0.6276,   0.6276]]


#### CH3Cl
##a = [
##["C" ,  0     ,   0     ,   0     ], 
##["Cl",  0     ,   0     ,   1.7810],
##["H" ,  1.0424,   0     ,  -0.3901],
##["H" , -0.5212,   0.9027,  -0.3901],
##["H" , -0.5212,  -0.9027,  -0.3901]]


#### C6H6
##a = [
##["C",    0.00000000,   1.39499067,   0.00000000],
##["C",   -1.20809736,   0.69749534,   0.00000000],
##["C",   -1.20809736,  -0.69749534,   0.00000000],
##["C",    0.00000000,  -1.39499067,   0.00000000],
##["C",    1.20809736,  -0.69749533,   0.00000000],
##["C",    1.20809736,   0.69749534,   0.00000000],
##["H",    0.00000000,   2.49460097,   0.00000000],
##["H",   -2.16038781,   1.24730049,   0.00000000],
##["H",   -2.16038781,  -1.24730049,   0.00000000],
##["H",    0.00000000,  -2.49460097,   0.00000000],
##["H",    2.16038781,  -1.24730048,   0.00000000],
##["H",    2.16038781,   1.24730049,   0.00000000]]
###############################################################################


if __name__ == "__main__":

    a = [
         ["C", -1.2131, -0.6884, 0.0000],
         ["C", -1.2028, 0.7064, 0.0001],
         ["C", -0.0103, -1.3948, 0.0000],
         ["C", 0.0104, 1.3948, -0.0001],
         ["C", 1.2028, -0.7063, 0.0000],
         ["C", 1.2131, 0.6884, 0.0000],
         ["H", -2.1577, -1.2244, 0.0000],
         ["H", -2.1393, 1.2564, 0.0001],
         ["H", -0.0184, -2.4809, -0.0001],
         ["H", 0.0184, 2.4808, 0.0000],
         ["H", 2.1394, -1.2563, 0.0001],
         ["H", 2.1577, 1.2245, 0.0000]
    ]

    I, D = main(a)
    print("The inertia_tensor I is \n{:}\n".format(I))
    print("The diagonalized I is \n{:}\n".format(D))