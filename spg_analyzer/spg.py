from tabnanny import check
import data
from sym_op import *
from utils import *
import numpy as np
from data import INERTIA_TOLERANCE, TOLERANCE, DEG_TOLERANCE, DEGENERACY_TOLERANCE

#! Define class



# Molecule class, which has the following attributes:
# 1. mol: a 2D list of atom names and coordinates (cartesian) (like pyscf)
# atom name: string, coordinates: array of floats
# Example:
#  [['H', [0.0, 0.0, 0.0]], 
#   ['H', [0.0, 0.0, 1.0]],
#   ['O', [0.0 1.0 0.0]]]

# 2. natm: number of atoms in the molecule e.g. 3
# 3. atm_name: list of atom names  ['H', 'H', 'O'] 
# 4. atm_num: atomic numbers sequence  [1, 1, 8]
# 5. coordinates: list of atom coordinates [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]



class Molecule:
    '''
    Molecule class, which has the following attributes:
    1. name: name of the molecule
    2. natm: number of atoms in the molecule
    3. atm_name: list of atom names
    4. atm_num: atomic numbers sequence
    5. atm_mass: atomic mass sequence
    6. coordinates: list of atom coordinates
    7. headblock: header block of the molecule
    8. atmblock: atom block of the molecule
    9. bondblocktoend: bond block to the end of the molecule

    '''
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
                raise ValueError(f"Atom {atm} not valid.")
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
    '''
    SPG class, inherited from Molecule class, which has the following attributes:
    1. mol: Molecule object
    2. inertia: inertia tensor of the molecule
    3. prin_mmt: principal moments of inertia
    4. prin_axes: principal axes of inertia
    5. degeneracy: degeneracy of the molecule
    6. sym_type: symmetry type of the molecule
    '''

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
        if self.mol.natm == 1:
            self.spg = 'Kh'
            return

        self.inertia_tensor()

        # principal moments and axes of inertia
        self.principal_moments_axes()

        # align the principal axes of inertia with the x, y, z axes
        self.align_axes()

        # check the degeneracy of the molecule
        self.check_degeneracy()
        # print("degeneracy", self.degeneracy)

        # check the symmetry type of the molecule
        self.check_sym_type()
        # print("symmetry type", self.sym_type)

        # check the symmetry of the molecule
        if self.sym_type == "linear":
            self.check_symmetry_linear()
        elif self.sym_type == "asymmetric":
            self.check_symmetry_asymmetric()
        elif self.sym_type == "symmetric":
            self.check_symmetry_symmetric()
        elif self.sym_type == "spherical":
            self.check_symmetry_spherical()

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
        self.inertia = tensor

        return 
    
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

        self.prin_mmt = thre_cut(eig_val)
        print("principal moment", self.prin_mmt)
        self.prin_axes = thre_cut(eig_vec)
        print("principal axes", self.prin_axes)

        return

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
                rot_mat[i][j] = np.matmul(self.prin_axes[i], xyz_axis[j])
        
        # align the principal axes with the x, y, z axes
        self.mol.coordinates = np.matmul(self.mol.coordinates, rot_mat)
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
            self.degeneracy = 3
        elif np.count_nonzero(mmt_diff > tol) == 2:
            self.degeneracy = 2
        elif np.count_nonzero(mmt_diff > tol) == 3:
            self.degeneracy = 1
        else:
            self.degeneracy = 1

        return

    def check_sym_type(self, tol=TOLERANCE):
        '''
        Check the symmetry type of the molecule. (linear, symmetric, asymmetric, spherical)

        Linear:
        IA = IB != 0, IC = 0

        Symmetric (degeneracy = 2):
        IA = IB != IC

        Asymmetric (degeneracy = 1):
        IA != IB != IC

        Spherical (degeneracy = 3):
        IA = IB = IC

        Returns:
            type: string

        '''
        # [IA-IB, IB-IC, IA-IC]
        mmt = sorted(self.prin_mmt)

        if abs(mmt[0]) < tol and mmt[1] == mmt[2]:
            print(self.mol.name, "is a linear molecule.")
            self.sym_type =  "linear"
        if self.degeneracy == 3:
            print(self.mol.name, "is a spherical molecule.")
            self.sym_type =  "spherical"
        if self.degeneracy == 2:
            print(self.mol.name, "is a symmetric molecule.")
            self.sym_type =  "symmetric"
        if self.degeneracy == 1:
            print(self.mol.name, "is an asymmetric molecule.")
            self.sym_type =  "asymmetric"
        return

    def check_symmetry_linear(self):
        '''
        Check the symmetry of the linear molecule.

        If the molecule has a inversion center, it is Dooh, otherwise it is Coov.

        Symmetric elements:
        Dooh: E 2C∞ ∞sigma_i i 2S∞ ∞C2
        Coov: E 2C∞ sigma_v
        '''

        inv = check_inversion(self.mol.coordinates)
        if inv:
            return "Dooh"
        else:
            return "Coov"

    def check_symmetry_asymmetric(self):
        '''
        Check the symmetry of the asymmetric molecule.
        '''



        return

    def check_symmetry_symmetric(self):
        '''
        Check the symmetry of the symmetric molecule.
        '''



        return

    def check_symmetry_spherical(self):
        '''
        Check the symmetry of the spherical molecule.

        Possible point groups:
        T: T, Th, Td
        O: O, Oh
        I: I, Ih

        Symmetric elements:
        T: E 4C3 4C3 3C2
        Th: E 4C3 4C3 3C2 i 4S6 4S6 3sigma_h
        Td: E 8C3 3C2 6S4 6sigma_d
        O: E 6C4 3C2 8C3 6C2
        Oh: E 8C3 6C2 6C4 3C2 i 6S4 8S6 3sigma_h 6sigma_d
        I: E 12C5 12C5 20C3 15C2
        Ih: E 12C5 12C5 20C3 15C2 i 12S10 12S10 20S6 15sigma

            T  Th  Td  O  Oh  I  Ih
        i   n  y   n   n   y   n   y
        Cn  3  3   3   4   4   5   5
        Sn  n  6   4   n   6   n   10

        '''


        # First check inversion
        # Yes -> Ih, Oh, Th, No -> I, O, Td, Th
        # Check the highest Cn
        # 3 -> T, 4 -> O, 5 -> I
        # Check the highest Sn
        # 6 -> Th, 4 -> Td

        Cn = max(check_Cn_extd(self.mol.coordinates), check_Cn(self.mol.coordinates))
        # print("Cn", Cn)

        if check_inversion(self.mol.coordinates):
            if Cn == 5:
                self.spg = "Ih"
            elif Cn == 4:
                self.spg = "Oh"
            elif Cn == 3:
                self.spg = "Th"
        else:
            if Cn == 5:
                self.spg = "I"
            elif Cn == 4:
                self.spg = "O"
            elif Cn == 3:
                Sn = max(check_Sn_extd(self.mol.coordinates), check_Sn(self.mol.coordinates))
                # print("Sn", Sn)
                if Sn == 6:
                    self.spg = "Th"
                elif Sn == 4:
                    self.spg = "Td"
                else:
                    self.spg = "T"
            else:
                self.spg = "TBD"

        return

# if __name__ == "__main__":
#     test_mol()