
from spg import Molecule


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
                if inv_atm[j][0:3]==atm[k].tolist() and inv_atm[j][3]==atm[k+1]:
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