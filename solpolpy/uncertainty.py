import copy
import numpy as np
import astropy.units as u

from solpolpy.alpha import radial_north
from ndcube import NDCollection, NDCube

def BpBuncertainty(input_collection):
    # Using the theoretical Error estimations, taking pB and pB' for ideal MZP
    alpha = input_collection["alpha"].data * u.radian
    pB, pBp = input_collection["pB"].data, input_collection["pBp"].data
    
    dth = (MZP_ang - off_MZP) *np.pi*u.radian/(180*u.degree)
    
    dpB = 4/3 * pBp * np.sqrt((np.cos(4*(MZP_ang[0]-alpha))*dth[0])**2 + (np.cos(4*(MZP_ang[1]-alpha))*dth[1])**2 + (np.cos(4*(MZP_ang[2]-alpha))*dth[2])**2)
    dpBp = 4/3 * pB * np.sqrt((np.cos(4*(MZP_ang[0]-alpha))*dth[0])**2 + (np.cos(4*(MZP_ang[1]-alpha))*dth[1])**2 + (np.cos(4*(MZP_ang[2]-alpha))*dth[2])**2)
    dB = np.sqrt(2/3)*(pB**2 + pBp**2)*dth[0]

    metaB, metapB, metapBp = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta), copy.copy(
                                input_collection["pBp"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB'), metapBp.update(Polar='pB-prime')
    Bp3_uncert_cube = []
    Bp3_uncert_cube.append(("B_uncert", NDCube(dB, wcs=input_collection["B"].wcs,  meta=metaB)))
    Bp3_uncert_cube.append(("pB_uncert", NDCube(dpB, wcs=input_collection["B"].wcs,  meta=metapB)))
    Bp3_uncert_cube.append(("pBp_uncert", NDCube(dpBp, wcs=input_collection["B"].wcs,  meta=metapBp)))
    Bp3_uncert_cube.append(("alpha", NDCube(alpha, wcs=input_collection["B"].wcs)))

    return NDCollection(Bp3_uncert_cube, meta={}, aligned_axes="all")