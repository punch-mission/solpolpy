import numpy as np
import copy
import astropy.units as u
from ndcube import NDCube, NDCollection


def conv_polar_from_head(input_cube):
    return int(float(str(input_cube.meta['POLAR']).strip(" Deg")))

#TODO: prepare a config file where the reference angle say of STEREO, KCor etc can be set
def npol_to_mzp(input_cube):
    """
    Notes
    ------
    Equation 44 in DeForest et al. 2022.

    """""
    input_dict = {}
    in_list = list(input_cube)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    # constants come from https://www.sciencedirect.com/science/article/pii/S0019103515003620?via%3Dihub
    if input_cube['angle_1'].meta['OBSRVTRY'] == 'STEREO_B':
        offset_angle = -18 * u.degree # STEREOB
    elif input_cube['angle_1'].meta['OBSRVTRY'] == 'STEREO_A':
        offset_angle = 45.8 * u.degree # STEREOA
    else:
        offset_angle = 0

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(conv_polar_from_head(input_cube[p_angle])) * u.degree * conv_fact] = input_cube[p_angle].data

    mzp_ang = [-60, 0, 60]
    Bmzp = {}
    for ang in mzp_ang: Bmzp[ang * u.degree] = (1 / 3) * np.sum(
        [ith_polarizer_brightness * (1 + 2 * np.cos(2 * (ang * u.degree * conv_fact - (ith_angle-offset_angle))))
         for ith_angle, ith_polarizer_brightness in input_dict.items()], axis=0)

    # todo: update header properly; time info?
    metaM, metaZ, metaP = copy.copy(input_cube["angle_1"].meta), copy.copy(input_cube["angle_2"].meta), copy.copy(input_cube["angle_3"].meta)
    metaM.update(POLAR=-60), metaZ.update(POLAR=0), metaP.update(POLAR=60)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_cube["angle_1"].wcs, meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_cube["angle_1"].wcs, meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_cube["angle_1"].wcs, meta=metaP)))
    for p_angle in in_list:
        if p_angle.lower() == "alpha":
            Bmzp_cube.append(("alpha", NDCube(input_cube['alpha'].data * u.radian, wcs=input_cube["angle_1"].wcs, meta=metaP)))
    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def mzp_to_bpb(input_cube):
    """
    Notes
    ------
    Equation 7 and 9 in DeForest et al. 2022.

    """""
    # TODO: need to check if 3 angles are input.
    # TODO: need to check if separated appropriately if not create quality warning.
    input_dict = {}
    in_list = list(input_cube)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_cube[p_angle].meta['POLAR']) * u.degree * conv_fact] = input_cube[p_angle].data

    alpha = input_cube['alpha'].data * u.radian
    B = (2 / 3) * (np.sum([ith_polarizer_brightness
                           for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))
    # todo: update header properly
    metaB, metapB = copy.copy(input_cube["Bm"].meta), copy.copy(input_cube["Bm"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB')

    BpB_cube = []
    BpB_cube.append(("B", NDCube(B, wcs=input_cube["Bm"].wcs, meta=metaB)))
    BpB_cube.append(("pB", NDCube(pB, wcs=input_cube["Bm"].wcs, meta=metapB)))
    BpB_cube.append(("alpha", NDCube(alpha, wcs=input_cube["Bm"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


def bpb_to_mzp(input_cube):
    """
    Notes
    ------
    Equation 4 in DeForest et al. 2022.
    """

    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    alpha = input_cube['alpha'].data * u.radian
    B, pB = input_cube["B"].data, input_cube["pB"].data
    mzp_ang = [-60, 0, 60]
    Bmzp = {}
    for ang in mzp_ang: Bmzp[ang * u.degree] = (1 / 2) * (B - pB * (np.cos(2 * (ang * u.degree - alpha))))

    metaM, metaZ, metaP = copy.copy(input_cube["B"].meta), copy.copy(input_cube["B"].meta), copy.copy(
        input_cube["B"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_cube["B"].wcs, meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_cube["B"].wcs, meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_cube["B"].wcs, meta=metaP)))
    Bmzp_cube.append(("alpha", NDCube(alpha, wcs=input_cube["B"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def bpb_to_btbr(input_cube):
    """
    Notes
    ------
    Equation 1 and 2 in DeForest et al. 2022.
    """
    input_dict = {}
    in_list = list(input_cube)

    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_cube[p_angle].meta['POLAR'])] = input_cube[p_angle].data

    alpha = input_cube['alpha'].data * u.radian
    B, pB = input_dict['B'], input_dict['pB']
    Br = (B - pB) / 2
    Bt = (B + pB) / 2

    metaBr, metaBt = copy.copy(input_cube["B"].meta), copy.copy(input_cube["B"].meta)
    metaBr.update(Polar='Br'), metaBt.update(Polar='Bt')

    BtBr_cube = []
    BtBr_cube.append(("Bt", NDCube(Bt, wcs=input_cube["B"].wcs, meta=metaBt)))
    BtBr_cube.append(("Br", NDCube(Br, wcs=input_cube["B"].wcs, meta=metaBr)))
    BtBr_cube.append(("alpha", NDCube(alpha, wcs=input_cube["B"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BtBr_cube, meta={}, aligned_axes="all")


def btbr_to_bpb(input_cube):
    """
    Notes
    ------
    Equation in Table 1 in DeForest et al. 2022.
    """
    input_dict = {}
    in_list = list(input_cube)

    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_cube[p_angle].meta['POLAR'])] = input_cube[p_angle].data

    alpha = input_cube['alpha'].data * u.radian
    Bt, Br = input_dict['Bt'], input_dict['Br']
    pB = (Bt - Br)
    B = (Bt + Br)

    metaB, metapB = copy.copy(input_cube["Bt"].meta), copy.copy(input_cube["Bt"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB')

    BpB_cube = []
    BpB_cube.append(("B", NDCube(B, wcs=input_cube["Bt"].wcs, meta=metaB)))
    BpB_cube.append(("pB", NDCube(pB, wcs=input_cube["Bt"].wcs, meta=metapB)))
    BpB_cube.append(("alpha", NDCube(alpha, wcs=input_cube["Bt"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


def mzp_to_stokes(input_cube):
    """
    Notes
    ------
    Equation 9, 12 and 13 in DeForest et al. 2022.
    """
    Bm, Bz, Bp = input_cube["Bm"].data, input_cube["Bz"].data, input_cube["Bp"].data

    mulmx = (2 / 3) * np.array([[1, 1, 1], [-1, 2, -1], [-np.sqrt(3), 0, np.sqrt(3)]])

    Bi = mulmx[0, 0] * Bm + mulmx[0, 1] * Bz + mulmx[0, 2] * Bp
    Bq = mulmx[1, 0] * Bm + mulmx[1, 1] * Bz + mulmx[1, 2] * Bp
    Bu = mulmx[2, 0] * Bm + mulmx[2, 1] * Bz + mulmx[2, 2] * Bp

    metaI, metaQ, metaU = copy.copy(input_cube["Bm"].meta), copy.copy(input_cube["Bz"].meta), copy.copy(
        input_cube["Bp"].meta)
    metaI.update(Polar='Stokes I'), metaQ.update(Polar='Stokes Q'), metaU.update(Polar='Stokes U')
    BStokes_cube = []
    BStokes_cube.append(("Bi", NDCube(Bi, wcs=input_cube["Bm"].wcs, meta=metaI)))
    BStokes_cube.append(("Bq", NDCube(Bq, wcs=input_cube["Bm"].wcs, meta=metaQ)))
    BStokes_cube.append(("Bu", NDCube(Bu, wcs=input_cube["Bm"].wcs, meta=metaU)))
    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")


def stokes_to_mzp(input_cube):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022. with alpha = np.pi/2
    """

    alpha = np.pi / 2
    Bi, Bq, Bu = input_cube["Bi"].data, input_cube["Bq"].data, input_cube["Bu"].data

    inv_mul_mx = (1 / 2) * np.array([[1, -np.cos(2 * (-np.pi / 3 - alpha)), -np.sin(2 * (-np.pi / 3 - alpha))],
                                     [1, -np.cos(2 * (0 - alpha)), 0],
                                     [1, -np.cos(2 * (np.pi / 3 - alpha)), -np.sin(2 * (np.pi / 3 - alpha))]])

    Bm = inv_mul_mx[0, 0] * Bi + inv_mul_mx[0, 1] * Bq + inv_mul_mx[0, 2] * Bu
    Bz = inv_mul_mx[1, 0] * Bi + inv_mul_mx[1, 1] * Bq + inv_mul_mx[1, 2] * Bu
    Bp = inv_mul_mx[2, 0] * Bi + inv_mul_mx[2, 1] * Bq + inv_mul_mx[2, 2] * Bu

    metaM, metaZ, metaP = copy.copy(input_cube["Bi"].meta), copy.copy(input_cube["Bq"].meta), copy.copy(
        input_cube["Bu"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bm, wcs=input_cube["Bi"].wcs, meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bz, wcs=input_cube["Bi"].wcs, meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bp, wcs=input_cube["Bi"].wcs, meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def mzp_to_bp3(input_cube):
    """
    Notes
    ------
    Equation 7, 9 and 10 in DeForest et al. 2022.
    """""
    input_dict = {}
    in_list = list(input_cube)
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[(input_cube[p_angle].meta['POLAR'] * u.degree * conv_fact)] = input_cube[p_angle].data
        
    alpha = input_cube['alpha'].data * u.radian
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
    metaB, metapB, metapBp = copy.copy(input_cube["Bm"].meta), copy.copy(input_cube["Bm"].meta), copy.copy(
        input_cube["Bm"].meta)
    metaB.update(Polar='B'), metapB.update(Polar='pB'), metapBp.update(Polar='pB-prime')

    Bp3_cube = []
    Bp3_cube.append(("B", NDCube(B, wcs=input_cube["Bm"].wcs, meta=metaB)))
    Bp3_cube.append(("pB", NDCube(pB, wcs=input_cube["Bm"].wcs, meta=metapB)))
    Bp3_cube.append(("pBp", NDCube(pBp, wcs=input_cube["Bm"].wcs, meta=metapBp)))
    Bp3_cube.append(("alpha", NDCube(alpha, wcs=input_cube["Bm"].wcs)))
    # todo: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bp3_cube, meta={}, aligned_axes="all")


def bp3_to_mzp(input_cube):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022.
    """""
    conv_fact = (np.pi * u.radian) / (180 * u.degree)

    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    B, pB, pBp = input_cube['B'].data, input_cube['pB'].data, input_cube['pBp'].data
    alpha = input_cube['alpha'].data * u.radian

    mzp_ang = [-60, 0, 60] * u.degree
    Bmzp = {}
    for ang in mzp_ang: Bmzp[ang] = (1 / 2) * (B - np.cos(2 * (ang * conv_fact - alpha)) * pB -
                                               np.cos(2 * (ang * conv_fact - alpha)) * pBp)

    metaM, metaZ, metaP = copy.copy(input_cube["B"].meta), copy.copy(input_cube["pB"].meta), copy.copy(
        input_cube["pBp"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_cube["B"].wcs, meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_cube["B"].wcs, meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_cube["B"].wcs, meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def btbr_to_mzp(input_cube):
    """
    Notes
    ------
    Equation 3 in DeForest et al. 2022.
    """
    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    alpha = input_cube['alpha'].data * u.radian
    Bt, Br = input_cube['Bt'].data, input_cube['Br'].data

    mzp_ang = [-60, 0, 60] * u.degree
    Bmzp = {}
    for ang in mzp_ang: Bmzp[ang] = Bt * (np.sin(ang - alpha)) ** 2 + Br * (np.cos(ang - alpha)) ** 2

    metaM, metaZ, metaP = copy.copy(input_cube["Bt"].meta), copy.copy(input_cube["Bt"].meta), copy.copy(
        input_cube["Bt"].meta)
    metaM.update(Polar=-60), metaZ.update(Polar=0), metaP.update(Polar=60)
    Bmzp_cube = []
    Bmzp_cube.append(("Bm", NDCube(Bmzp[-60 * u.degree], wcs=input_cube["Bt"].wcs, meta=metaM)))
    Bmzp_cube.append(("Bz", NDCube(Bmzp[0 * u.degree], wcs=input_cube["Bt"].wcs, meta=metaZ)))
    Bmzp_cube.append(("Bp", NDCube(Bmzp[60 * u.degree], wcs=input_cube["Bt"].wcs, meta=metaP)))

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


def bp3_to_bthp(input_cube):
    """
    Notes
    ------
    Equations 9, 15, 16 in DeForest et al. 2022.
    """""
    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    B, pB, pBp = input_cube['B'].data, input_cube['pB'].data, input_cube['pBp'].data
    alpha = input_cube['alpha'].data * u.radian

    theta_mx = (1 / 2) * np.arctan2(pBp, pB) * u.radian + np.pi / 2 * u.radian + alpha
    p = np.sqrt(pB ** 2 + pBp ** 2) / B

    metaTh, metaP = copy.copy(input_cube["B"].meta), copy.copy(input_cube["pB"].meta)
    metaTh.update(Polar='Theta'), metaP.update(Polar='Degree of Polarization')

    Bthp_cube = []
    Bthp_cube.append(("B", NDCube(B, wcs=input_cube["B"].wcs, meta=input_cube['B'].meta)))
    Bthp_cube.append(("theta", NDCube(theta_mx, wcs=input_cube["B"].wcs, meta=metaTh)))
    Bthp_cube.append(("p", NDCube(p, wcs=input_cube["B"].wcs, meta=metaP)))

    return NDCollection(Bthp_cube, meta={}, aligned_axes="all")


def btbr_to_npol(input_cube, angles):
    """
    Notes
    ------
    Equation 3 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    if "alpha" not in input_cube:
        raise ValueError("missing alpha")

    alpha = input_cube['alpha'].data * u.radian
    Bt, Br = input_cube['Bt'].data, input_cube['Br'].data

    npol_ang = angles
    Bnpol = {}
    Bnpol_cube = []
    for ang in npol_ang:
        Bnpol[ang] = Bt * (np.sin(ang * u.degree - alpha)) ** 2 + Br * (np.cos(ang * u.degree - alpha)) ** 2
        meta_tmp = copy.copy(input_cube["Bt"].meta)
        meta_tmp.update(Polar=(ang))
        Bnpol_cube.append(('B' + str(ang), NDCube(Bnpol[ang], wcs=input_cube["Bt"].wcs, meta=meta_tmp)))
    Bnpol_cube.append(("alpha", NDCube(alpha, wcs=input_cube["Bt"].wcs)))

    return NDCollection(Bnpol_cube, meta={}, aligned_axes="all")


def fourpol_to_stokes(input_cube):
    """
    Notes
    ------
    Table 1 in DeForest et al. 2022.

    """""
    Bi = input_cube["B0"].data + input_cube["B90"].data
    Bq = input_cube["B90"].data - input_cube["B0"].data
    Bu = input_cube["B135"].data - input_cube["B45"].data

    metaI, metaQ, metaU = copy.copy(input_cube["B0"].meta), copy.copy(input_cube["B0"].meta), copy.copy(input_cube["B0"].meta)
    BStokes_cube = []
    BStokes_cube.append(("Bi", NDCube(Bi, wcs=input_cube["B0"].wcs, meta=metaI)))
    BStokes_cube.append(("Bq", NDCube(Bq, wcs=input_cube["B0"].wcs, meta=metaQ)))
    BStokes_cube.append(("Bu", NDCube(Bu, wcs=input_cube["B0"].wcs, meta=metaU)))

    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")
