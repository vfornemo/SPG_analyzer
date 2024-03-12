from sympy import deg
import mol_parser, data, sym_op
from utils import *
import numpy as np
from data import INERTIA_TOLERANCE, TOLERANCE, DEG_TOLERANCE, DEGENERACY_TOLERANCE

#! Define class



# Define a class called Molecule, which has the following attributes:
# 1. mol: a 2D list of atom names and coordinates (cartesian) (like pyscf)
# atom name: string, coordinates: array of floats
# Example:
#  [['H' [0.0, 0.0, 0.0]], 
#   ['H' [0.0, 0.0, 1.0]],
#   ['O' [0.0 1.0 0.0]]]

# 2. natm: number of atoms in the molecule e.g. 3
# 3. atm_name: list of atom names  ['H', 'H', 'O'] 
# 4. atm_num: atomic numbers sequence  [1, 1, 8]
# 5. coordinates: list of atom coordinates [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]



class Molecule:
    def __init__(self, mol):
        self.name = "mol"
        # number of atoms
        self.natm = len(mol)
        # atomic name sequence
        self.atm_name = [atom[0] for atom in mol]
        # coordinates
        self.coordinates = [atom[1] for atom in mol]

    def build(self):
        '''
        Build the molecule.
        '''
        # atomic number sequence
        self.atm_num = self.get_atm_num()
        # atomic mass sequence
        self.atm_mass = self.get_atm_mass()

        # move the molecule to the mass center
        mass_ctr = self.mass_center()
        self.coord_shift(mass_ctr)
        self.coordinates = thre_cut(self.coordinates)
        # print("coordinates", self.coordinates)

        return

    
    # Get atomic number from atomic symbol
    def get_atm_num(self):
        num_list = []
        for atm in self.atm_name:
            if atm in data.atom_dict:
                num_list.append(data.atom_dict[atm])
            else:
                RuntimeError(f"Atom {atm} not valid.")
        return num_list
    
    # Get atomic mass from atomic number
    def get_atm_mass(self):
        mass_list = []
        for atm in self.atm_num:
            mass_list.append(data.atom_data[atm][3])
        return mass_list

    # Get mass center of the molecule
    def mass_center(self):
        mass = 0
        center = [0, 0, 0]
        for i in range(self.natm):
            mass += self.atm_mass[i]
            for j in range(3):
                center[j] += self.atm_mass[i] * self.coordinates[i][j]
        for j in range(3):
            center[j] /= mass
        return center
    
    # Shift the coordinates of the molecule
    def coord_shift(self, shift):
        for i in range(self.natm):
            for j in range(3):
                self.coordinates[i][j] -= shift[j]
        return

    # def atm_block_updater(self):
    #     '''
    #     Update the atom block of the molecule.
    #     '''
    #     for i in range(self.natm):
    #         for j in range(3):
    #             self.atmblock[i][j] = format(self.coordinates[i][j], '.4f')
    #     return



# Define a class called SPG(symmetric point group), inherited from Molecule class, which has the following attributes:
# This class is used to analyze the point group of a molecule
# 1. name: name of the point group
# 2. mol: molecule object
    

class SPG(Molecule):
    def __init__(self, mol):
        '''
        Args:
            mol: Molecule object
        '''
        self.mol = mol

    def build(self):
        '''
        Build the data for symmetry analysis, including computing the inertia tensor, principal moments and axes of inertia, 
        and aligning the principal axes of inertia with the x, y, z axes.
        '''
        self.inertia = self.inertia_tensor()
        # principal moments and axes of inertia
        self.prin_mmt, self.prin_axes = self.principal_moments_axes()
        self.prin_mmt = thre_cut(self.prin_mmt)
        print("principal moment", self.prin_mmt)
        self.prin_axes = thre_cut(self.prin_axes)
        print("principal axes", self.prin_axes)

        # align the principal axes of inertia with the x, y, z axes
        self.align_axes()

        # check the degeneracy of the molecule
        self.degeneracy = self.check_degeneracy()

        # check the symmetry type of the molecule
        self.sym_type = self.check_sym_type()



        return

    def inertia_tensor(self, tol = INERTIA_TOLERANCE):
        '''
        Calculate the inertia tensor of the molecule.
        [[Ixx, Ixy, Ixz],
         [Iyx, Iyy, Iyz],
         [Izx, Izy, Izz]]   
        '''
        tensor = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(self.mol.natm):
            for j in range(3):
                for k in range(3):
                    # Sigma(m_i * (delta(j,k) * r_i^2 - r_i[j] * r_i[k]))
                    tensor[j][k] += self.mol.atm_mass[i] * (delta(j,k)*distance_square(self.mol.coordinates[i]) - self.mol.coordinates[i][j] * self.mol.coordinates[i][k])
        tol = int(-np.log10(tol))
        tensor = thre_cut(tensor, tol)
        print("inertia tensor after threshold cut", tensor)
        return tensor
    
    def principal_moments_axes(self):
        '''
        Calculate the principal moments and axes of inertia of the molecule.
        moments: [I1 I2 I3]
        axes: [[n1x n1y n1z], x
               [n2x n2y n2z], y
               [n3x n3y n3z]] z
        '''
    	# Diagonalize the inertia tensor
        eig_val, eig_vec = np.linalg.eig(self.inertia)
        # eig_val is the principal moments of inertia
        # eig_vec is the principal axes of inertia
        
        return eig_val, eig_vec	

    def align_axes(self, tol = DEG_TOLERANCE):
        '''
        Align the principal axes of inertia with the x, y, z axes.
        '''
        # rotation matrix
        # [[cosxx', cosxy', cosxz']
        #  [cosyx', cosyy', cosyz']
        #  [coszx', coszy', coszz']]

        # set the tolerance for the alignment, angle between the principal axis and x, y, z axes
        tol = np.cos(np.deg2rad(tol))
        # print("tol", tol)

        # calculate the angle theta between the principal axis and z axis
        xyz_axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        # rotation matrix
        # X = MX', M' = M^T
        # X' = M^T X
        rot_mat = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                rot_mat[i][j] = np.dot(self.prin_axes[i], xyz_axis[j])
        
        # align the principal axes with the x, y, z axes
        self.mol.coordinates = np.dot(self.mol.coordinates, rot_mat)
        self.mol.coordinates = thre_cut(self.mol.coordinates)
        return
    
    def check_degeneracy(self, tol=DEGENERACY_TOLERANCE):
        '''
        Check the degeneracy of the molecule.
        
        Symmetric (degeneracy = 2):
        IA = IB != IC
        Asymmetric (degeneracy = 1):
        IA != IB != IC
        Spherical (degeneracy = 3):
        IA = IB = IC

        '''
        # [IA-IB, IB-IC, IA-IC]
        mmt_diff = np.array([abs(self.prin_mmt[0] - self.prin_mmt[1]), abs(self.prin_mmt[1] - self.prin_mmt[2]), abs(self.prin_mmt[0] - self.prin_mmt[2])])
        if np.count_nonzero(mmt_diff > tol) == 0:
            return 3
        elif np.count_nonzero(mmt_diff > tol) == 2:
            return 2
        elif np.count_nonzero(mmt_diff > tol) == 3:
            return 1
        else:
            return 1

    def check_sym_type(self, tol=TOLERANCE):
        '''
        Check the symmetry type of the molecule. (linear, symmetric, asymmetric, spherical)

        Linear:
        IA = IB != 0, IC = 0

        Returns:
            type: string

        '''
        # [IA-IB, IB-IC, IA-IC]
        mmt = sorted(self.prin_mmt)

        if abs(mmt[0]) < tol and mmt[1] == mmt[2]:
            print(self.mol.name, "is a linear molecule.")
            return "linear"
        if self.degeneracy == 3:
            print(self.mol.name, "is a spherical molecule.")
            return "spherical"
        if self.degeneracy == 2:
            print(self.mol.name, "is a symmetric molecule.")
            return "symmetric"
        if self.degeneracy == 1:
            print(self.mol.name, "is an asymmetric molecule.")
            return "asymmetric"
        return
        



    def check_symmetry(self):
        '''
        Check the symmetry of the molecule.
        '''

        if self.mol.natm == 1:
            self.spg = 'R3'
            return
        
        # if self.prin_mmt


        return


# test to check if the function works

def test_mol():
    water = mol_parser.from_mol('../tests/testset/4_Cs.mol')
    # water = mol_parser.from_mol('../tests/mol_samples/H2O2.mol')

    water.build()
    water_spg = SPG(water)
    mol_parser.to_mol(water, water.name + '_original')
    water_spg.build()
    mol_parser.to_mol(water_spg.mol, water.name + '_aligned')

    return

if __name__ == "__main__":
    test_mol()