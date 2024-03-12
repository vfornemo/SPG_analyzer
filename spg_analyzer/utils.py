"""
This module, `utils.py`, provides a set of utility functions for performing calculations related to points in 3D space. 

Functions:
- `delta(i, j)`: Implements the Kronecker delta.
- `distance(coord1, coord2)`: Calculates the Euclidean distance between two points in 3D space.
- `distance_square(coord1, coord2)`: Similar to `distance`, but returns the square of the distance.
- `car2sph(coords)`: Converts a list of points from Cartesian coordinates to spherical coordinates.
- `sph2car(coords)`: Converts a list of points from spherical coordinates to Cartesian coordinates.

"""

import numpy as np

from data import TOLERANCE

# Kronecker delta
def delta(i, j):
    if i == j:
        return 1
    else:
        return 0

def distance(coord1, coord2=[0, 0, 0]):
    '''
    Calculate the distance between two points.
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

def thre_cut(x, thre=6):
    '''
    Cut the number to given decimal places.

    Args:
        x: array
        thre: threshold for rounding, number of decimal places.
    
    '''
    return np.around(x, thre)
