"""
Polarizer functions for solpolpy.
Found in:
DeForest, C. E., Seaton, D. B., & West, M. J. (2022).
Three-polarizer Treatment of Linear Polarization in Coronagraphs and Heliospheric Imagers.
The Astrophysical Journal, 927(1), 98.
"""
import copy

import astropy.units as u
import numpy as np
from ndcube import NDCollection, NDCube


def conv_polar_from_head(input_cube):
    return int(float(str(input_cube.meta['POLAR']).strip(" Deg")))


# TODO: prepare a config file where the reference angle say of STEREO, KCor etc can be set
def npol_to_mzp(input_collection):
    """
    Notes
    ------
    Equation 44 in DeForest et al. 2022.

    """""
    input_dict = {}
    in_list = list(input_collection)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    # constants come from https://www.sciencedirect.com/science/article/pii/S0019103515003620?via%3Dihub
    if input_collection['angle_1'].meta['OBSRVTRY'] == 'STEREO_B':
        offset_angle = -18 * u.degree  # STEREOB
    elif input_collection['angle_1'].meta['OBSRVTRY'] == 'STEREO_A':
        offset_angle = 45.8 * u.degree  # STEREOA
    else:
        offset_angle = 0

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(conv_polar_from_head(input_collection[p_angle])) * u.degree * conv_fact] = input_collection[p_angle].data

    mzp_ang = [-60, 0, 60]
    Bmzp = {}
    for ang in mzp_ang:
        Bmzp[ang * u.degree] = ((1 / 3)
                                * np.sum([ith_polarizer_brightness
                                          * (1 + 2 * np.cos(2 * (ang * u.degree * conv_fact
                                                                 - (ith_angle-offset_angle))))
                                          for ith_angle, ith_polarizer_brightness in input_dict.items()], axis=0))

    # todo: update header properly; time info?
    metaM, metaZ, metaP = (copy.copy(input_collection["angle_1"].meta),
                           copy.copy(input_collection["angle_2"].meta),
                           copy.copy(input_collection["angle_3"].meta))
    metaM.update(POLAR=-60), metaZ.update(POLAR=0), metaP.update(POLAR=60)
    mask = combine_masks(input_collection["angle_1"].mask,
                                 input_collection["angle_2"].mask,
                                 input_collection["angle_3"].mask)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["angle_1"].wcs, mask=mask,  meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_collection["angle_1"].wcs, mask=mask,  meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_collection["angle_1"].wcs, mask=mask,  meta=metaP)))
    for p_angle in in_list:
        if p_angle.lower() == "alpha":
            Bmzp_cube.append(("alpha", NDCube(input_collection['alpha'].data * u.radian,
                                              wcs=input_collection["angle_1"].wcs,
                                              meta=metaP)))
    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def mzp_to_bpb(input_collection):
    """
    Notes
    ------
    Equation 7 and 9 in DeForest et al. 2022.

    """""
    # TODO: need to check if 3 angles are input.
    # TODO: need to check if separated appropriately if not create quality warning.
    input_dict = {}
    in_list = list(input_collection)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_collection[p_angle].meta['POLAR']) * u.degree * conv_fact] = input_collection[p_angle].data

    alpha = input_collection['alpha'].data * u.radian
    B = (2 / 3) * (np.sum([ith_polarizer_brightness
                           for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))
    metaB, metapB = copy.copy(input_collection["Bm"].meta), copy.copy(input_collection["Bm"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB')
    mask = combine_masks(input_collection["Bm"].mask, input_collection["Bz"].mask, input_collection["Bp"].mask)
    BpB_cube = []
    BpB_cube.append(("B", NDCube(B, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metaB)))
    BpB_cube.append(("pB", NDCube(pB, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metapB)))
    BpB_cube.append(("alpha", NDCube(alpha, wcs=input_collection["Bm"].wcs, mask=mask)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


def bpb_to_mzp(input_collection):
    """
    Notes
    ------
    Equation 4 in DeForest et al. 2022.
    """

    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    alpha = input_collection['alpha'].data * u.radian
    B, pB = input_collection["B"].data, input_collection["pB"].data
    mzp_ang = [-60, 0, 60]
    Bmzp = {}
    for ang in mzp_ang:
        Bmzp[ang * u.degree] = (1 / 2) * (B - pB * (np.cos(2 * (ang * u.degree - alpha))))

    metaM, metaZ, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["B"].meta), copy.copy(
        input_collection["B"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    mask = combine_masks(input_collection["B"].mask, input_collection["pB"].mask)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaP)))
    Bmzp_cube.append(("alpha", NDCube(alpha, wcs=input_collection["B"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def bpb_to_btbr(input_collection):
    """
    Notes
    ------
    Equation 1 and 2 in DeForest et al. 2022.
    """

    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    alpha = input_collection['alpha'].data * u.radian
    B, pB = input_collection['B'].data, input_collection['pB'].data
    Br = (B - pB) / 2
    Bt = (B + pB) / 2

    metaBr, metaBt = copy.copy(input_collection["B"].meta), copy.copy(input_collection["B"].meta)
    metaBr.update(Polar='Br'), metaBt.update(Polar='Bt')
    mask = combine_masks(input_collection["B"].mask, input_collection["pB"].mask)
    BtBr_cube = []
    BtBr_cube.append(("Bt", NDCube(Bt, wcs=input_collection["B"].wcs, mask=mask,  meta=metaBt)))
    BtBr_cube.append(("Br", NDCube(Br, wcs=input_collection["B"].wcs, mask=mask,  meta=metaBr)))
    BtBr_cube.append(("alpha", NDCube(alpha, wcs=input_collection["B"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BtBr_cube, meta={}, aligned_axes="all")


def btbr_to_bpb(input_collection):
    """
    Notes
    ------
    Equation in Table 1 in DeForest et al. 2022.
    """

    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    alpha = input_collection['alpha'].data * u.radian
    Bt, Br = input_collection['Bt'].data, input_collection['Br'].data
    pB = (Bt - Br)
    B = (Bt + Br)

    metaB, metapB = copy.copy(input_collection["Bt"].meta), copy.copy(input_collection["Bt"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB')
    mask = combine_masks(input_collection["Bt"].mask, input_collection["Br"].mask)
    BpB_cube = []
    BpB_cube.append(("B", NDCube(B, wcs=input_collection["Bt"].wcs, mask=mask,  meta=metaB)))
    BpB_cube.append(("pB", NDCube(pB, wcs=input_collection["Bt"].wcs, mask=mask,  meta=metapB)))
    BpB_cube.append(("alpha", NDCube(alpha, wcs=input_collection["Bt"].wcs, mask=mask)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


def mzp_to_stokes(input_collection):
    """
    Notes
    ------
    Equation 9, 12 and 13 in DeForest et al. 2022.
    """
    Bm, Bz, Bp = input_collection["Bm"].data, input_collection["Bz"].data, input_collection["Bp"].data

    mulmx = (2 / 3) * np.array([[1, 1, 1], [-1, 2, -1], [-np.sqrt(3), 0, np.sqrt(3)]])

    Bi = mulmx[0, 0] * Bm + mulmx[0, 1] * Bz + mulmx[0, 2] * Bp
    Bq = mulmx[1, 0] * Bm + mulmx[1, 1] * Bz + mulmx[1, 2] * Bp
    Bu = mulmx[2, 0] * Bm + mulmx[2, 1] * Bz + mulmx[2, 2] * Bp

    metaI, metaQ, metaU = copy.copy(input_collection["Bm"].meta), copy.copy(input_collection["Bz"].meta), copy.copy(
        input_collection["Bp"].meta)
    metaI.update(Polar='Stokes I'), metaQ.update(Polar='Stokes Q'), metaU.update(Polar='Stokes U')
    mask = combine_masks(input_collection["Bm"].mask, input_collection["Bz"].mask, input_collection["Bp"].mask)
    BStokes_cube = []
    BStokes_cube.append(("Bi", NDCube(Bi, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metaI)))
    BStokes_cube.append(("Bq", NDCube(Bq, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metaQ)))
    BStokes_cube.append(("Bu", NDCube(Bu, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metaU)))
    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")


def stokes_to_mzp(input_collection):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022. with alpha = np.pi/2
    """

    alpha = np.pi / 2
    Bi, Bq, Bu = input_collection["Bi"].data, input_collection["Bq"].data, input_collection["Bu"].data

    inv_mul_mx = (1 / 2) * np.array([[1, -np.cos(2 * (-np.pi / 3 - alpha)), -np.sin(2 * (-np.pi / 3 - alpha))],
                                     [1, -np.cos(2 * (0 - alpha)), 0],
                                     [1, -np.cos(2 * (np.pi / 3 - alpha)), -np.sin(2 * (np.pi / 3 - alpha))]])

    Bm = inv_mul_mx[0, 0] * Bi + inv_mul_mx[0, 1] * Bq + inv_mul_mx[0, 2] * Bu
    Bz = inv_mul_mx[1, 0] * Bi + inv_mul_mx[1, 1] * Bq + inv_mul_mx[1, 2] * Bu
    Bp = inv_mul_mx[2, 0] * Bi + inv_mul_mx[2, 1] * Bq + inv_mul_mx[2, 2] * Bu

    metaM, metaZ, metaP = copy.copy(input_collection["Bi"].meta), copy.copy(input_collection["Bq"].meta), copy.copy(
        input_collection["Bu"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    mask = combine_masks(input_collection["Bi"].mask, input_collection["Bq"].mask, input_collection["Bu"].mask)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bm, wcs=input_collection["Bi"].wcs, mask=mask,  meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bz, wcs=input_collection["Bi"].wcs, mask=mask,  meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bp, wcs=input_collection["Bi"].wcs, mask=mask,  meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def mzp_to_bp3(input_collection):
    """
    Notes
    ------
    Equation 7, 9 and 10 in DeForest et al. 2022.
    """""
    input_dict = {}
    in_list = list(input_collection)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_collection[p_angle].meta['POLAR'] * u.degree * conv_fact)] = input_collection[p_angle].data

    alpha = input_collection['alpha'].data * u.radian
    B = (2 / 3) * (np.sum([ith_polarizer_brightness for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))

    pBp = (-4 / 3) * (np.sum([ith_polarizer_brightness * np.sin(2 * (ith_angle - alpha))
                              for ith_angle, ith_polarizer_brightness
                              in input_dict.items() if ith_angle != "alpha"], axis=0))
    # todo: update header properly
    metaB, metapB, metapBp = copy.copy(input_collection["Bm"].meta), copy.copy(input_collection["Bm"].meta), copy.copy(
        input_collection["Bm"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB'), metapBp.update(Polar='pB-prime')
    mask = combine_masks(input_collection["Bm"].mask, input_collection["Bz"].mask, input_collection["Bp"].mask)
    Bp3_cube = []
    Bp3_cube.append(("B", NDCube(B, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metaB)))
    Bp3_cube.append(("pB", NDCube(pB, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metapB)))
    Bp3_cube.append(("pBp", NDCube(pBp, wcs=input_collection["Bm"].wcs, mask=mask,  meta=metapBp)))
    Bp3_cube.append(("alpha", NDCube(alpha, wcs=input_collection["Bm"].wcs,  mask=mask)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bp3_cube, meta={}, aligned_axes="all")


def bp3_to_mzp(input_collection):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022.
    """""
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    B, pB, pBp = input_collection['B'].data, input_collection['pB'].data, input_collection['pBp'].data
    alpha = input_collection['alpha'].data * u.radian

    mzp_ang = [-60, 0, 60] * u.degree
    Bmzp = {}
    for ang in mzp_ang:
        Bmzp[ang] = (1 / 2) * (B - np.cos(2 * (ang * conv_fact - alpha)) * pB -
                               np.cos(2 * (ang * conv_fact - alpha)) * pBp)

    metaM, metaZ, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta), copy.copy(
        input_collection["pBp"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    mask = combine_masks(input_collection["B"].mask, input_collection["pB"].mask, input_collection["pBp"].mask)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_collection["B"].wcs, mask=mask,  meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def btbr_to_mzp(input_collection):
    """
    Notes
    ------
    Equation 3 in DeForest et al. 2022.
    """
    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    alpha = input_collection['alpha'].data * u.radian
    Bt = input_collection['Bt'].data
    Br = input_collection['Br'].data

    mzp_ang = [-60, 0, 60] * u.degree
    Bmzp = {}
    for ang in mzp_ang:
        Bmzp[ang] = Bt * (np.sin(ang - alpha)) ** 2 + Br * (np.cos(ang - alpha)) ** 2

    metaM, metaZ, metaP = copy.copy(input_collection["Bt"].meta), copy.copy(input_collection["Bt"].meta), copy.copy(
        input_collection["Bt"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    mask = combine_masks(input_collection["Bt"].mask, input_collection["Br"].mask)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask,  meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask,  meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask,  meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def bp3_to_bthp(input_collection):
    """
    Notes
    ------
    Equations 9, 15, 16 in DeForest et al. 2022.
    """""
    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    B, pB, pBp = input_collection['B'].data, input_collection['pB'].data, input_collection['pBp'].data
    alpha = input_collection['alpha'].data * u.radian

    theta_mx = (1 / 2) * np.arctan2(pBp, pB) * u.radian + np.pi / 2 * u.radian + alpha
    p = np.sqrt(pB ** 2 + pBp ** 2) / B

    metaTh, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta)
    metaTh.update(Polar='Theta'), metaP.update(Polar='Degree of Polarization')
    mask = combine_masks(input_collection["B"].mask, input_collection["pB"].mask, input_collection["pBp"].mask)
    Bthp_cube = []
    Bthp_cube.append(("B", NDCube(B, wcs=input_collection["B"].wcs, mask=mask,  meta=input_collection['B'].meta)))
    Bthp_cube.append(("theta", NDCube(theta_mx, wcs=input_collection["B"].wcs, mask=mask,  meta=metaTh)))
    Bthp_cube.append(("p", NDCube(p, wcs=input_collection["B"].wcs, mask=mask,  meta=metaP)))

    return NDCollection(Bthp_cube, meta={}, aligned_axes="all")


def btbr_to_npol(input_collection, angles):
    """
    Notes
    ------
    Equation 3 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    if "alpha" not in input_collection:
        raise ValueError("missing alpha")

    alpha = input_collection['alpha'].data * u.radian
    Bt, Br = input_collection['Bt'].data, input_collection['Br'].data

    npol_ang = angles
    Bnpol = {}
    Bnpol_cube = []
    mask = combine_masks(input_collection["Bt"].mask, input_collection["Br"].mask)
    for ang in npol_ang:
        Bnpol[ang] = Bt * (np.sin(ang * u.degree - alpha)) ** 2 + Br * (np.cos(ang * u.degree - alpha)) ** 2
        meta_tmp = copy.copy(input_collection["Bt"].meta)
        meta_tmp.update(Polar=(ang))
        Bnpol_cube.append(('B' + str(ang), NDCube(Bnpol[ang], wcs=input_collection["Bt"].wcs, mask=mask,  meta=meta_tmp)))
    Bnpol_cube.append(("alpha", NDCube(alpha, wcs=input_collection["Bt"].wcs, mask=mask)))

    return NDCollection(Bnpol_cube, meta={}, aligned_axes="all")


def fourpol_to_stokes(input_collection):
    """
    Notes
    ------
    Table 1 in DeForest et al. 2022.

    """""
    Bi = input_collection["B0"].data + input_collection["B90"].data
    Bq = input_collection["B90"].data - input_collection["B0"].data
    Bu = input_collection["B135"].data - input_collection["B45"].data

    metaI, metaQ, metaU = (copy.copy(input_collection["B0"].meta),
                           copy.copy(input_collection["B0"].meta),
                           copy.copy(input_collection["B0"].meta))
    mask = combine_masks(input_collection["B0"].mask,
                                 input_collection["B45"].mask,
                                 input_collection["B90"].mask,
                                 input_collection["B135"].mask)
    BStokes_cube = []
    BStokes_cube.append(("Bi", NDCube(Bi, wcs=input_collection["B0"].wcs, mask=mask,  meta=metaI)))
    BStokes_cube.append(("Bq", NDCube(Bq, wcs=input_collection["B0"].wcs, mask=mask,  meta=metaQ)))
    BStokes_cube.append(("Bu", NDCube(Bu, wcs=input_collection["B0"].wcs, mask=mask,  meta=metaU)))
    BStokes_cube.append(("Bu", NDCube(Bu, wcs=input_collection["B0"].wcs, mask=mask,  meta=metaU)))

    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")


def combine_masks(*args):
    """ Combines masks

    If any of the masks are None, the result is None.
    Otherwise, when combining any value that is masked in any of the input args, gets masked, i.e. it does a logical or.
    """
    if any(arg is None for arg in args):
        return None
    else:
        return np.logical_or.reduce(args)
