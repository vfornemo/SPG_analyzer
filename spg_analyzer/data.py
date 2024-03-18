"""
This file contains the element data for the spg_analyzer package.

"""

from termios import VLNEXT
from zmq import MECHANISM


TOLERANCE = 1e-5
DEG_TOLERANCE = 4.0
INERTIA_TOLERANCE = 1e-3
DEGENERACY_TOLERANCE = 1e-2

VERY_TIGHT_TOLERANCE = 1e-4
TIGHT_TOLERANCE = 1e-3
MEDIUM_TOLERANCE = 1e-2
LOOSE_TOLERANCE = 1e-1
VERY_LOOSE_TOLERANCE = 0.3
SUPER_LOOSE_TOLERANCE = 0.5
ULTRA_LOOSE_TOLERANCE = 1.0

tol_map = {
    "very_tight": VERY_TIGHT_TOLERANCE,
    "tight": TIGHT_TOLERANCE,
    "medium": MEDIUM_TOLERANCE,
    "loose": LOOSE_TOLERANCE,
    "very_loose": VERY_LOOSE_TOLERANCE,
    "super_loose": SUPER_LOOSE_TOLERANCE,
    "ultra_loose": ULTRA_LOOSE_TOLERANCE,
}


atom_data = [
    # atomic number, symbols, names, masses, bohr radius
    [  0, "X", "X",            0.000000, 0.000],  # 0
    [  1, "H", "Hydrogen",     1.007940, 0.324],  # 1
    [  2, "He", "Helium",      4.002602, 0.000],  # 2
    [  3, "Li", "Lithium",     6.941000, 1.271],  # 3
    [  4, "Be", "Beryllium",   9.012182, 0.927],  # 4
    [  5, "B", "Boron",       10.811000, 0.874],  # 5
    [  6, "C", "Carbon",      12.010700, 0.759],  # 6
    [  7, "N", "Nitrogen",    14.006700, 0.706],  # 7
    [  8, "O", "Oxygen",      15.999400, 0.678],  # 8
    [  9, "F", "Fluorine",    18.998403, 0.568],  # 9
    [ 10, "Ne", "Neon",       20.179700, 0.000],  # 10
    [ 11, "Na", "Sodium",     22.989769, 1.672],  # 11
    [ 12, "Mg", "Magnesium",  24.305000, 1.358],  # 12
    [ 13, "Al", "Aluminum",  26.981539, 1.218],  # 13
    [ 14, "Si", "Silicon",    28.085500, 1.187],  # 14
    [ 15, "P", "Phosphorus",  30.973762, 1.105],  # 15
    [ 16, "S", "Sulfur",      32.065000, 1.045],  # 16
    [ 17, "Cl", "Chlorine",   35.453000, 1.006],  # 17
    [ 18, "Ar", "Argon",      39.948000, 0.000],  # 18
    [ 19, "K", "Potassium",   39.098300, 2.247],  # 19
    [ 20, "Ca", "Calcium",    40.078000, 1.748],  # 20
    [ 21, "Sc", "Scandium",   44.955912, 1.664],  # 21
    [ 22, "Ti", "Titanium",   47.867000, 1.620],  # 22
    [ 23, "V", "Vanadium",    50.941500, 1.543],  # 23
    [ 24, "Cr", "Chromium",   51.996100, 1.418],  # 24
    [ 25, "Mn", "Manganese",  54.938045, 1.569],  # 25
    [ 26, "Fe", "Iron",       55.845000, 1.514],  # 26
    [ 27, "Co", "Cobalt",     58.933195, 1.385],  # 27
    [ 28, "Ni", "Nickel",     58.693400, 1.390],  # 28
    [ 29, "Cu", "Copper",     63.546000, 1.382],  # 29
    [ 30, "Zn", "Zinc",       65.380000, 1.416],  # 30
    [ 31, "Ga", "Gallium",    69.723000, 1.235],  # 31
    [ 32, "Ge", "Germanium",  72.640000, 1.201],  # 32
    [ 33, "As", "Arsenic",    74.921600, 1.232],  # 33
    [ 34, "Se", "Selenium",   78.960000, 1.210],  # 34
    [ 35, "Br", "Bromine",    79.904000, 1.190],  # 35
    [ 36, "Kr", "Krypton",    83.798000, 0.000],  # 36
    [ 37, "Rb", "Rubidium",   85.467800, 2.284],  # 37
    [ 38, "Sr", "Strontium",  87.620000, 1.942],  # 38
    [ 39, "Y", "Yttrium",     88.905850, 1.993],  # 39
    [ 40, "Zr", "Zirconium",  91.224000, 1.758],  # 40
    [ 41, "Nb", "Niobium",    92.906380, 1.610],  # 41
    [ 42, "Mo", "Molybdenum", 95.960000, 1.639],  # 42
    [ 43, "Tc", "Technetium",  0.000000, 1.493],  # 43
    [ 44, "Ru", "Ruthenium",  101.07000, 1.467],  # 44
    [ 45, "Rh", "Rhodium",    102.90550, 1.437],  # 45
    [ 46, "Pd", "Palladium",  106.42000, 1.422],  # 46
    [ 47, "Ag", "Silver",     107.86820, 1.466],  # 47
    [ 48, "Cd", "Cadmium",    112.41100, 1.441],  # 48
    [ 49, "In", "Indium",     114.81800, 1.421],  # 49
    [ 50, "Sn", "Tin",        118.71000, 1.408],  # 50
    [ 51, "Sb", "Antimony",   121.76000, 1.397],  # 51
    [ 52, "Te", "Tellurium",  127.60000, 1.395],  # 52
    [ 53, "I", "Iodine",      126.90447, 1.396],  # 53
    [ 54, "Xe", "Xenon",      131.29300, 1.336],  # 54
    [ 55, "Cs", "Caesium",    132.90545, 2.470],  # 55
    [ 56, "Ba", "Barium",     137.32700, 2.219],  # 56
    [ 57, "La", "Lanthanum",  138.90547, 2.089],  # 57
    [ 58, "Ce", "Cerium",     140.11600, 2.054],  # 58
    [ 59, "Pr", "Praseodymium", 140.90765, 1.979],  # 59
    [ 60, "Nd", "Neodymium",  144.24200, 0.000],  # 60
    [ 61, "Pm", "Promethium",   0.00000, 0.000],  # 61
    [ 62, "Sm", "Samarium",   150.36000, 2.535],  # 62
    [ 63, "Eu", "Europium",   151.96400, 0.000],  # 63
    [ 64, "Gd", "Gadolinium", 157.25000, 0.000],  # 64
    [ 65, "Tb", "Terbium",    158.92535, 0.000],  # 65
    [ 66, "Dy", "Dysprosium", 162.50000, 0.000],  # 66
    [ 67, "Ho", "Holmium",    164.93032, 0.000],  # 67
    [ 68, "Er", "Erbium",     167.25900, 0.000],  # 68
    [ 69, "Tm", "Thulium",    168.93421, 0.000],  # 69
    [ 70, "Yb", "Ytterbium",  173.05400, 0.000],  # 70
    [ 71, "Lu", "Lutetium",   174.96680, 0.000],  # 71
    [ 72, "Hf", "Hafnium",    178.49000, 1.779],  # 72
    [ 73, "Ta", "Tantalum",   180.94788, 1.723],  # 73
    [ 74, "W", "Tungsten",    183.84000, 1.627],  # 74
    [ 75, "Re", "Rhenium",    186.20700, 1.536],  # 75
    [ 76, "Os", "Osmium",     190.23000, 1.521],  # 76
    [ 77, "Ir", "Iridium",    192.21700, 1.456],  # 77
    [ 78, "Pt", "Platinum",   195.08400, 1.390],  # 78
    [ 79, "Au", "Gold",       196.96657, 1.402],  # 79
    [ 80, "Hg", "Mercury",    200.59000, 1.371],  # 80
    [ 81, "Tl", "Thallium",   204.38330, 1.384],  # 81
    [ 82, "Pb", "Lead",       207.20000, 1.820],  # 82
    [ 83, "Bi", "Bismuth",    208.98040, 1.507],  # 83
    [ 84, "Po", "Polonium",     0.00000, 0.000],  # 84
    [ 85, "At", "Astatine",     0.00000, 0.000],  # 85
    [ 86, "Rn", "Radon",        0.00000, 0.000],  # 86
    [ 87, "Fr", "Francium",     0.00000, 0.000],  # 87
    [ 88, "Ra", "Radium",       0.00000, 0.000],  # 88
    [ 89, "Ac", "Actinium",     0.00000, 0.000],  # 89
    [ 90, "Th", "Thorium",    232.03806, 0.000],  # 90
    [ 91, "Pa", "Protactinium",231.03588, 0.000], # 91
    [ 92, "U", "Uranium",     238.02891, 0.000],  # 92
    [ 93, "Np", "Neptunium",    0.00000, 0.000],  # 93
    [ 94, "Pu", "Plutonium",    0.00000, 0.000],  # 94
    [ 95, "Am", "Americium",    0.00000, 0.000],  # 95
    [ 96, "Cm", "Curium",       0.00000, 0.000],  # 96
    [ 97, "Bk", "Berkelium",    0.00000, 0.000],  # 97
    [ 98, "Cf", "Californium",  0.00000, 0.000],  # 98
    [ 99, "Es", "Einsteinium",  0.00000, 0.000],  # 99
    [100, "Fm", "Fermium",      0.00000, 0.000],  # 100
    [101, "Md", "Mendelevium",  0.00000, 0.000],  # 101
    [102, "No", "Nobelium",     0.00000, 0.000],  # 102
    [103, "Lr", "Lawrencium",   0.00000, 0.000],  # 103
    [104, "Rf", "Rutherfordium",0.00000, 0.000],  # 104
    [105, "Db", "Dubnium",      0.00000, 0.000],  # 105
    [106, "Sg", "Seaborgium",   0.00000, 0.000],  # 106
    [107, "Bh", "Bohrium",      0.00000, 0.000],  # 107
    [108, "Hs", "Hassium",      0.00000, 0.000],  # 108
    [109, "Mt", "Meitnerium",   0.00000, 0.000],  # 109
    [110, "Ds", "Darmstadtium", 0.00000, 0.000],  # 110
    [111, "Rg", "Roentgenium",  0.00000, 0.000],  # 111
    [112, "Cn", "Copernicium",  0.00000, 0.000],  # 112
    [113, "Nh", "Nihonium",   0.00000, 0.000],  # 113
    [114, "Fl", "Flerovium", 0.00000, 0.000],  # 114
    [115, "Mc", "Moscovium", 0.00000, 0.000],  # 115
    [116, "Lv", "Livermorium",  0.00000, 0.000],  # 116
    [117, "Ts", "Tennessine", 0.00000, 0.000],  # 117
    [118, "Og", "Oganesson",  0.00000, 0.000],  # 118
    ]

# Hashmap of atomic symbol to atomic number 1-118
atom_dict = {
    'X'  : 0,
    'H'  : 1,
    'He' : 2,
    'Li' : 3,
    'Be' : 4,
    'B'  : 5,
    'C'  : 6,
    'N'  : 7,
    'O'  : 8,
    'F'  : 9,
    'Ne' : 10,
    'Na' : 11,
    'Mg' : 12,
    'Al' : 13, 
    'Si' : 14,
    'P'  : 15,
    'S'  : 16,
    'Cl' : 17,
    'Ar' : 18,
    'K'  : 19,
    'Ca' : 20,
    'Sc' : 21,
    'Ti' : 22,
    'V'  : 23,
    'Cr' : 24,
    'Mn' : 25,
    'Fe' : 26,
    'Co' : 27,
    'Ni' : 28,
    'Cu' : 29,
    'Zn' : 30,
    'Ga' : 31,
    'Ge' : 32,
    'As' : 33,
    'Se' : 34,
    'Br' : 35,
    'Kr' : 36,
    'Rb' : 37,
    'Sr' : 38,
    'Y'  : 39,
    'Zr' : 40,
    'Nb' : 41,
    'Mo' : 42,
    'Tc' : 43,
    'Ru' : 44,
    'Rh' : 45,
    'Pd' : 46,
    'Ag' : 47,
    'Cd' : 48,
    'In' : 49,
    'Sn' : 50,
    'Sb' : 51,
    'Te' : 52,
    'I'  : 53,
    'Xe' : 54,
    'Cs' : 55,
    'Ba' : 56,
    'La' : 57,
    'Ce' : 58,
    'Pr' : 59,
    'Nd' : 60,
    'Pm' : 61,
    'Sm' : 62,
    'Eu' : 63,
    'Gd' : 64,
    'Tb' : 65,
    'Dy' : 66,
    'Ho' : 67,
    'Er' : 68,
    'Tm' : 69,
    'Yb' : 70,
    'Lu' : 71,
    'Hf' : 72,
    'Ta' : 73,
    'W'  : 74,
    'Re' : 75,
    'Os' : 76,
    'Ir' : 77,
    'Pt' : 78,
    'Au' : 79,
    'Hg' : 80,
    'Tl' : 81,
    'Pb' : 82,
    'Bi' : 83,
    'Po' : 84,
    'At' : 85,
    'Rn' : 86,
    'Fr' : 87,
    'Ra' : 88,
    'Ac' : 89,
    'Th' : 90,
    'Pa' : 91,
    'U'  : 92,
    'Np' : 93,
    'Pu' : 94,
    'Am' : 95,
    'Cm' : 96,
    'Bk' : 97,
    'Cf' : 98,
    'Es' : 99,
    'Fm' : 100,
    'Md' : 101,
    'No' : 102,
    'Lr' : 103,
    'Rf' : 104,
    'Ha' : 105,
    'Sg' : 106,
    'Bh' : 107,
    'Hs' : 108,
    'Mt' : 109,
    'Ds' : 110,
    'Rg' : 111,
    'Cn' : 112,
    'Nh' : 113,
    'Fl' : 114,
    'Mc' : 115,
    'Lv' : 116,
    'Ts' : 117,
    'Og' : 118,
}

# Extended axis for finding the symmetry operations
# Including 1 the axis between the x, y, and z axix
#          [[ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
#           [-1,  1, 0], [-1, 0,  1], [0, -1,  1],
#           [ 1, -1, 0], [ 1, 0, -1], [0,  1, -1],
#           [-1, -1, 0], [-1, 0, -1], [0, -1, -1],
#           2 the axis among the x, y, and z axix
#           [1, 1, 1], 
#           [-1, 1, 1], [1, -1, 1], [1, 1, -1], 
#           [-1, -1, 1], [-1, 1, -1], [1, -1, -1], 
#           [-1, -1, -1]]
# Actually, [1, 1, 1] is the same as [-1, -1, -1], so we can simplify the list

extend_axis = [[ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
                     [-1,  1, 0], [-1, 0,  1], [0, -1,  1],
                     [1, 1, 1], 
                     [-1, 1, 1], [1, -1, 1], [1, 1, -1]]

full_axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [ 1,  1, 0], [ 1, 0,  1], [0,  1,  1], 
                     [-1,  1, 0], [-1, 0,  1], [0, -1,  1],
                     [1, 1, 1], 
                     [-1, 1, 1], [1, -1, 1], [1, 1, -1]]

