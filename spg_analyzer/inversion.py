
from spg import *

at1=[0.8110,-0.4677,-0.0000,"H"]
at2=[-0.0000,0.0589,0.0000,"O"]
at3=[-0.8110,-0.4677,-0.0000,"H"]


def invert(x, y, z, atom):
    xi=-x
    yi=-y
    zi=-z
    xyzi=[float("{:.4f}".format(xi)), float("{:.4f}".format(yi)), float("{:.4f}".format(zi)),"atom"]
    return xyzi

def test_inversion_H2O(c1,c2,c3):

    atom1=False
    atom2=False
    atom3=False
    i1=invert(at1[0],at1[1],at1[2],at1[3])
    i2=invert(at2[0],at2[1],at2[2],at2[3])
    i3=invert(at3[0],at3[1],at3[2],at3[3])

    if i1==at1 or i1==at2 or i1==at3:
        atom1=True
    if i2==at1 or i2==at2 or i2==at3:
        atom2=True
    if i3==at1 or i3==at2 or i3==at3:
        atom3=True
    
    if atom1==True and atom2==True and atom3==True:
        print("There is an inversion center")
        inv=True
    else:
        print("No inversion center")
        inv=False
    return inv
H2O=[['H' [0.0, 0.0, 0.0]],['H' [0.0, 0.0, 1.0]],['O' [0.0, 1.0, 0.0]]]


mol_test=Molecule(H2O,3,['H', 'H', 'O'],[1, 1, 8],[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
mol_test.natm


# def test_inversion(cd):

#     atom1=False
#     atom2=False
#     atom3=False
#     i1=invert(at1[0],at1[1],at1[2],at1[3])
#     i2=invert(at2[0],at2[1],at2[2],at2[3])
#     i3=invert(at3[0],at3[1],at3[2],at3[3])

#     if i1==at1 or i1==at2 or i1==at3:
#         atom1=True
#     if i2==at1 or i2==at2 or i2==at3:
#         atom2=True
#     if i3==at1 or i3==at2 or i3==at3:
#         atom3=True
    
#     if atom1==True and atom2==True and atom3==True:
#         print("There is an inversion center")
#         inv=True
#     else:
#         print("No inversion center")
#         inv=False
#     return inv


# inv1=test_inversion_H2O(at1,at2,at3)




    


