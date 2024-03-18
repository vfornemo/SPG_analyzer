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

class SPG(Molecule):
    '''
    SPG class, inherited from Molecule class, which has the following attributes:
    1. mol: Molecule object
    2. inertia: inertia tensor of the molecule
    3. prin_mmt: principal moments of inertia
    4. prin_axes: principal axes of inertia
    5. degeneracy: degeneracy of the molecule
    6. sym_type: symmetry type of the molecule
    7. so: symmetric operations of the molecule (distinguishable SOs for point group determination)
    8. spg: point group of the molecule
    9. mode: tolerance mode, 'ultra_loose', 'super_loose', 'loose', 'medium', 'tight', 'very_tight', default is 'medium'
    '''

    def __init__(self, mol):
        '''
        Args:
            mol: Molecule object
        '''
        self.mol = mol
        self.so = []
        self.mode = "medium"

    def build(self):
        '''
        Build the data for symmetry analysis, including computing the inertia tensor, principal moments and axes of inertia, 
        and aligning the principal axes of inertia with the x, y, z axes.
        '''
        if self.mol.natm == 1:
            self.spg = 'Kh'
            self.sym_type = 'spherical'
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
        self.check_symmetry()

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
        elif self.degeneracy == 3:
            print(self.mol.name, "is a spherical molecule.")
            self.sym_type =  "spherical"
        elif self.degeneracy == 2:
            print(self.mol.name, "is a symmetric molecule.")
            self.sym_type =  "symmetric"
        elif self.degeneracy == 1:
            print(self.mol.name, "is an asymmetric molecule.")
            self.sym_type =  "asymmetric"
        else:
            print("Symmetry type undetermined.")
            self.sym_type =  None
        return
    
    def check_symmetry(self):
        if self.sym_type == "linear":
            self.check_symmetry_linear()
        elif self.sym_type == "asymmetric":
            self.check_symmetry_asym()
        elif self.sym_type == "symmetric":
            self.check_symmetry_sym()
        elif self.sym_type == "spherical":
            self.check_symmetry_sph()
        else:
            print("Symmetry type undetermined.")
            self.spg = "TBD"
        return

    def check_symmetry_linear(self):
        '''
        Check the symmetry of the linear molecule.

        Possible point groups:
        Coov, Dooh

        If the molecule has a inversion center, it is Dooh, otherwise it is Coov.

        Symmetric elements:
        Dooh: E 2C∞ ∞sigma_i i 2S∞ ∞C2
        Coov: E 2C∞ sigma_v
        '''

        inv = check_inversion(self.mol.atm_name, self.mol.coordinates)
        if inv:
            self.so.append("i")
            self.spg = "Dooh"
        else:
            self.spg = "Coov"
        return

    def check_symmetry_asym(self):
        '''
        Check the symmetry of the asymmetric molecule.
        '''

        self.spg = "TBD"

        return

    def check_symmetry_sym(self):
        '''
        Check the symmetry of the symmetric molecule.

        Possible point groups:
        Dnh, Dnd, Dn, Cnh, Cnv, S2n

        Symmetric elements:

        S2n:

        Cnv:
        C3v: E 2C3 3sigma_v
        C4v: E 2C4 C2 2sigma_v 2sigma_d
        C5v: E 2C5 2C5 5sigma_v
        C6v: E 2C6 2C3 C2 3sigma_v 3sigma_d
        
        Cnh:
        C3h: E 2C3 sigma_h 2S3
        C4h: E 2C4 C2 i 2S4 sigma_h
        C5h: E 2C5 2C5 sigma_h 2S5 2S5
        C6h: E 2C6 2C3 C2 i 2S6 2S3 sigma_h

        Dn:
        D3: E 2C3 3C2
        D4: E 2C4 C2 2C2 2C2
        D5: E 2C5 2C5 5C2
        D6: E 2C6 2C3 C2 3C2 3C2

        Dnd:
        D3d: E 2C3 3C2 i 2S6 3sigma_d
        D4d: E 2S8 2C4 2S8 C2 4C2 4sigma_d
        D5d: E 2C5 2C5 5C2 i 2S10 2S10 5sigma_d
        D6d: E 2S12 2C6 2S4 2C3 2S12 C2 6C2 6sigma_d

        Dnh:
        D3h: E 2C3 3C2 sigma_h 2S3 3sigma_v
        D4h: E 2C4 C2 2C2 2C2 i 2S4 sigma_h 2sigma_v 2sigma_d
        D5h: E 2C5 2C5 5C2 sigma_h 2S5 2S5 5sigma_v
        D6h: E 2C6 2C3 C2 3C2 3C2 i 2S6 2S3 sigma_h 3sigma_v 3sigma_d
        
        '''
        Cn, self.Cn_axis = check_Cn(self.mol.atm_name, self.mol.coordinates, True, self.mode)
        nC2 = check_C2_perp_Cn(self.mol.atm_name, self.mol.coordinates, self.Cn_axis, self.mode)
        self.so.append("C" + str(Cn))
        if nC2:
            self.so.append(str(nC2) + "C2")

        if nC2 >= 1 and Cn >= 2:
            # Dnh Dnd Dn
            if check_reflection_h(self.mol.atm_name, self.mol.coordinates, self.Cn_axis, self.mode):
                # Dnh
                self.so.append("sigma_h")
                self.spg = "D" + str(Cn) + "h"
            else:
                # Dnd or Dn
                Sn = max(check_Sn_extd(self.mol.atm_name, self.mol.coordinates, 16, self.mode), check_Sn(self.mol.atm_name, self.mol.coordinates, 16, self.mode))
                if Sn:
                    self.so.append("S" + str(Sn))
                    self.spg = "D" + str(Cn) + "d"
                else:
                    self.spg = "D" + str(Cn)
        else:
            # Cnh Cnv S2n
            if check_reflection_h(self.mol.atm_name, self.mol.coordinates, self.Cn_axis, self.mode):
                # Cnh
                self.so.append("sigma_h")
                self.spg = "C" + str(Cn) + "h"
            else:
                # Cnv S2n
                Sn = max(check_Sn_extd(self.mol.atm_name, self.mol.coordinates, 8, self.mode), check_Sn(self.mol.atm_name, self.mol.coordinates, 8, self.mode))
                if Sn:
                    self.so.append("S" + str(Sn))
                    self.spg = "S" + str(2*Cn)
                else:
                    self.spg = "C" + str(Cn) + "v"

        return

    def check_symmetry_sph(self):
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

        Cn1, Cn_axis1 = check_Cn(self.mol.atm_name, self.mol.coordinates, True, self.mode)
        Cn2, Cn_axis2 = check_Cn_extd(self.mol.atm_name, self.mol.coordinates, True, self.mode)

        if Cn1 > Cn2:
            Cn = Cn1
            self.Cn_axis = Cn_axis1
        else:
            Cn = Cn2
            self.Cn_axis = Cn_axis2

        # print("Cn", Cn)

        if check_inversion(self.mol.atm_name, self.mol.coordinates, self.mode):
            self.so.append("i")
            if Cn == 5:
                self.so.append("C5")
                self.spg = "Ih"
            elif Cn == 4:
                self.so.append("C4")
                self.spg = "Oh"
            elif Cn == 3:
                self.so.append("C3")
                self.spg = "Th"
        else:
            if Cn == 5:
                self.so.append("C5")
                self.spg = "I"
            elif Cn == 4:
                self.so.append("C4")
                self.spg = "O"
            elif Cn == 3:
                self.so.append("C3")
                Sn = max(check_Sn_extd(self.mol.atm_name, self.mol.coordinates, 10, self.mode), check_Sn(self.mol.atm_name, self.mol.coordinates, 10, self.mode))
                # print("Sn", Sn)
                if Sn == 6:
                    self.so.append("S6")
                    self.spg = "TBD"
                elif Sn == 4:
                    self.so.append("S4")
                    self.spg = "Td"
                else:
                    self.spg = "T"
            else:
                self.spg = "TBD"

        return

# if __name__ == "__main__":
#     test_mol()