
from spg import Molecule


def invert(coord):
    xi=-coord[0]
    yi=-coord[1]
    zi=-coord[2]
    xyzi=[float("{:.4f}".format(xi)), float("{:.4f}".format(yi)), float("{:.4f}".format(zi)),coord[3]]
    return xyzi

#mol_test=Molecule(H2O,3,['H', 'H', 'O'],[1, 1, 8],[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
#mol1=Molecule('H' [0.0, 0.0, 0.0])

mol1 = Molecule([['H', [0.0, 0.0, 0.0]], 
   ['H', [0.0, 0.0, 1.0]],
   ['O', [0.0, 1.0, 0.0]]])

mol1.build()

# cd2=mol1.coordinates[0].tolist()
# cd2.append(mol1.atm_name[1])
# i2=invert(cd2,mol1.atm_name[2])
# old_atm1=cd2.append(mol1.atm_name[1])

# print(f"{cd2}")

# atm1=mol1.coordinates[0].tolist()
# atm1.append(mol1.atm_name[0])
# print(f"{atm1}")

def test_inversion(cd):
    atom1=False
    atom2=False
    atom3=False
    atm1=cd.coordinates[0].tolist()
    atm1.append(cd.atm_name[0])
    atm2=cd.coordinates[1].tolist()
    atm2.append(cd.atm_name[1])
    atm3=cd.coordinates[2].tolist()
    atm3.append(cd.atm_name[2])
    i1=invert(atm1)
    i2=invert(atm2)
    i3=invert(atm3)

    if i1==atm1 or i1==atm2 or i1==atm3:
        atom1=True
    if i2==atm1 or i2==atm2 or i2==atm3:
        atom2=True
    if i3==atm1 or i3==atm2 or i3==atm3:
        atom3=True
    
    if atom1==True and atom2==True and atom3==True:
        print("There is an inversion center")
        inv=True
    else:
        print("No inversion center")
        inv=False
    return inv



test_inversion(mol1)