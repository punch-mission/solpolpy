import copy

import astropy.units as u
import numpy as np
from ndcube import NDCollection, NDCube


# Polarizer misalignment
def polarizer_misalign(input_collection):
    # Using the theoretical Error estimations, taking pB and pB' for ideal MZP
    off_MZP = 0.1 * u.degree #TODO: Fix this using header information. Arbitrary value taken.
    alpha = input_collection["alpha"].data * u.radian
    pB, pBp = input_collection["pB"].data, input_collection["pBp"].data
    MZP_ang = [-60, 0, 60] * u.degree
    dth = (MZP_ang - off_MZP) * np.pi * u.radian / (180 * u.degree)

    dpB = 4 / 3 * pBp * np.sqrt(
        (np.cos(4 * (MZP_ang[0] - alpha)) * dth[0]) ** 2 + (np.cos(4 * (MZP_ang[1] - alpha)) * dth[1]) ** 2 + (
                np.cos(4 * (MZP_ang[2] - alpha)) * dth[2]) ** 2)
    dpBp = 4 / 3 * pB * np.sqrt(
        (np.cos(4 * (MZP_ang[0] - alpha)) * dth[0]) ** 2 + (np.cos(4 * (MZP_ang[1] - alpha)) * dth[1]) ** 2 + (
                np.cos(4 * (MZP_ang[2] - alpha)) * dth[2]) ** 2)
    dB = np.sqrt(2 / 3) * (pB ** 2 + pBp ** 2) * dth[0]

    metaB, metapB, metapBp = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta), copy.copy(
        input_collection["pBp"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB'), metapBp.update(Polar='pB-prime')
    Bp3_uncert_cube = [("B_uncert", NDCube(dB, wcs=input_collection["B"].wcs, meta=metaB)),
                       ("pB_uncert", NDCube(dpB, wcs=input_collection["B"].wcs, meta=metapB)),
                       ("pBp_uncert", NDCube(dpBp, wcs=input_collection["B"].wcs, meta=metapBp)),
                       ("alpha", NDCube(alpha, wcs=input_collection["B"].wcs))]

    return NDCollection(Bp3_uncert_cube, meta={}, aligned_axes="all")

# # Photometric error
# def photometric_uncert(input_collection):
#     alpha = input_collection["alpha"].data * u.radian
#     B, pB, pBp = input_collection["B"].data, input_collection["pB"].data, input_collection["pBp"].data
#
#     dpB = np.sqrt(dB)
#     dpBp = np.sqrt(dB)
#
#     metaB, metapB, metapBp = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta), copy.copy(
#         input_collection["pBp"].meta)
#     metaB.update(Polar='B'), metapB.update(Polar='pB'), metapBp.update(Polar='pB-prime')
#     Bp3_uncert_cube = [("B_uncert", NDCube(dB, wcs=input_collection["B"].wcs, meta=metaB)),
#                        ("pB_uncert", NDCube(dpB, wcs=input_collection["B"].wcs, meta=metapB)),
#                        ("pBp_uncert", NDCube(dpBp, wcs=input_collection["B"].wcs, meta=metapBp)),
#                        ("alpha", NDCube(alpha, wcs=input_collection["B"].wcs))]
#
#     return NDCollection(Bp3_uncert_cube, meta={}, aligned_axes="all")
