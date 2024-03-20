"""
This module, `utils.py`, provides a set of utility functions for performing calculations related to points in 3D space. 

Functions:
- `delta(i, j)`: Implements the Kronecker delta.
- `distance(coord1, coord2)`: Calculates the Euclidean distance between two points in 3D space.
- `distance_square(coord1, coord2)`: Similar to `distance`, but returns the square of the distance.
- `car2sph(coords)`: Converts a list of points from Cartesian coordinates to spherical coordinates.
- `sph2car(coords)`: Converts a list of points from spherical coordinates to Cartesian coordinates.

"""

from calendar import c
import numpy as np

from data import TOLERANCE, tol_map

# Kronecker delta
def delta(i, j):
    if i == j:
        return 1
    else:
        return 0

def distance(coord1, coord2=[0, 0, 0]):
    '''
    Calculate the distance between two points. If coord2 is not given, the distance from the origin is calculated.

    Args:
        coord1: list of [x, y, z]
        coord2: list of [x, y, z]   

    '''
    return ((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)**0.5

def distance_square(coord1, coord2=[0, 0, 0]):
    '''
    Calculate the distance between two points.
    '''
    return (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2

def car2sph(coords):
    '''
    Convert the cartesian coordinates to spherical coordinates.
    coords: list of [x, y, z]
    returns:
    sph_coords: list of [r, theta, phi]
    '''
    sph = np.zeros((len(coords),3))
    for i in range(len(coords)):
        sph[i][0] = (coords[i][0]**2 + coords[i][1]**2 + coords[i][2]**2)**0.5 # r
        sph[i][1] = np.arccos(coords[i][2]/sph[i][0]) # theta
        sph[i][2] = np.arctan2(coords[i][1], coords[i][0]) # phi
    return sph

def sph2car(coords):
    '''
    Convert the spherical coordinates to cartesian coordinates.
    coords: list of [r, theta, phi]
    returns:
    car_coords: list of [x, y, z]
    '''
    car = np.zeros((len(coords),3))
    for i in range(len(coords)):
        car[i][0] = coords[i][0]*np.sin(coords[i][1])*np.cos(coords[i][2]) # x
        car[i][1] = coords[i][0]*np.sin(coords[i][1])*np.sin(coords[i][2]) # y
        car[i][2] = coords[i][0]*np.cos(coords[i][1]) # z
    return car

def sort_coords_idx(coords):
    '''
    Sort the coordinates based on the x, y, z values.
    The sorting will first be based on the x values, then y values, then z values.
    This is useful for comparing the coordinates of the molecule before and after symmetric operations.
    '''
    sorted_idx = np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))

    return sorted_idx

def sort_coords(coords, sorted_idx):
    '''
    Return the sorted coordinates based on the sorted index.
    '''

    return coords[sorted_idx]

def sort_atoms(atm_list, sorted_idx):
    '''
    Return the sorted atom list based on the sorted coordinates.
    '''
  
    atm_list = np.array(atm_list)

    return atm_list[sorted_idx]

def thre_cut(x, thre=6):
    '''
    Cut the number to given decimal places.

    Args:
        x: array
        thre: threshold for rounding, number of decimal places.
    
    '''
    return np.around(x, thre)

def coords_intersect(atm_list1, atm_list2, coords1, coords2, mode = 'medium'):
    '''
    The sorted coordinates may have some twists, but they are equivalent. 
    This function checks if the two sets of coordinates are equivalent.

    Example:

    '''
    token = np.zeros(len(atm_list1))

    for i in range(len(atm_list1)):
        for j in range(len(atm_list2)):
            if np.allclose(coords1[i], coords2[j], atol=tol_map[mode]) and atm_list1[i] == atm_list2[j]:
                token[i] = 1
                break

    if np.all(token):
        return True
    else:
        return False 

# reference: https://en.wikipedia.org/wiki/3D_rotation_group, Mobius transformation
def rotate_molecule(coords):
    '''
    Rotate the molecule at random angles.
    
    Warning: rotating a molecule may reduce its symmetry, and the point group may change.
    Use GaussView or other software to check the symmetry of the molecule after rotation.

    Rotate matrix:
    R = [[cos(psi)*cos(phi)-cos(theta)*sin(psi)*sin(phi), -sin(psi)*cos(phi)-cos(theta)*cos(psi)*sin(phi), sin(theta)*sin(phi)],
         [cos(psi)*sin(phi)+cos(theta)*sin(psi)*cos(phi), -sin(psi)*sin(phi)+cos(theta)*cos(psi)*cos(phi), -sin(theta)*cos(phi)],
         [sin(theta)*sin(psi), sin(theta)*cos(psi), cos(theta)]]
    '''
    theta, phi, psi = np.random.rand(3)*2*np.pi
    R = np.array([[np.cos(psi)*np.cos(phi)-np.cos(theta)*np.sin(psi)*np.sin(phi), -np.sin(psi)*np.cos(phi)-np.cos(theta)*np.cos(psi)*np.sin(phi), np.sin(theta)*np.sin(phi)],
                    [np.cos(psi)*np.sin(phi)+np.cos(theta)*np.sin(psi)*np.cos(phi), -np.sin(psi)*np.sin(phi)+np.cos(theta)*np.cos(psi)*np.cos(phi), -np.sin(theta)*np.cos(phi)],
                    [np.sin(theta)*np.sin(psi), np.sin(theta)*np.cos(psi), np.cos(theta)]])
    
    return np.dot(R, coords.T).T



def compare_coords(atm_list, coords1, coords2, mode = 'medium', enhanced = False):
    '''
    Compare two sets of coordinates to see if they are the same within a tolerance.
    The coordinates are first sorted, then compared.

    Args:
        atm_list: list of atoms
        coords1: np.array of coordinates
        coords2: np.array of coordinates
        mode: string, the mode for tolerance, default is 'medium'
        enhanced: boolean, if True, the comparison will enable coords intersection function
                  to check if the two coords matrice are equivalent, default is False
    '''

    if mode == 'medium':
        cut = 3
    elif mode == 'tight':
        cut = 4
    elif mode == 'very_tight':
        cut = 5
    elif mode.__contains__('loose'):
        cut = 2
    else:
        cut = 3

    coords1 = thre_cut(coords1, cut)
    coords2 = thre_cut(coords2, cut)

    sorted_idx1 = sort_coords_idx(coords1)
    sorted_idx2 = sort_coords_idx(coords2)

    atm1 = sort_atoms(atm_list, sorted_idx1)
    atm2 = sort_atoms(atm_list, sorted_idx2)
    
    coords1 = sort_coords(coords1, sorted_idx1)
    coords2 = sort_coords(coords2, sorted_idx2)

    if enhanced:
        return coords_intersect(atm1, atm2, coords1, coords2, mode)
    else:
        return np.array_equal(atm1, atm2) and np.allclose(coords1, coords2, atol=tol_map[mode])


