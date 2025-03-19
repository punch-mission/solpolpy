"""Polarizer functions for solpolpy.
Found in:
DeForest, C. E., Seaton, D. B., & West, M. J. (2022).
Three-polarizer Treatment of Linear Polarization in Coronagraphs and Heliospheric Imagers.
The Astrophysical Journal, 927(1), 98.
"""
import sys
import copy
from enum import StrEnum
from inspect import signature, getmembers, isfunction

import astropy.units as u
import networkx as nx
import numpy as np
from ndcube import NDCollection, NDCube

from solpolpy.errors import InvalidDataError, MissingAlphaError, SolpolpyError
from solpolpy.util import combine_all_collection_masks, extract_crota_from_wcs

System = StrEnum("System", ["bpb", "npol", "stokes", "mzpsolar", "mzpinstru", "btbr", "bthp", "fourpol", "bp3"])
SYSTEM_REQUIRED_KEYS = {System.bpb: {"B", "pB"},
                        System.npol: set(),
                        System.stokes: {"I", "Q", "U"},
                        System.mzpsolar: {"M", "Z", "P"},
                        System.mzpinstru: {"M", "Z", "P"},
                        System.btbr: {"Bt", "Br"},
                        System.bp3: {"B", "pB", "pBp"},
                        System.bthp: {"B", "theta", "p"},
                        System.fourpol: {str(q) for q in [0.0, 45.0, 90.0, 135.0] * u.degree},
                        }


def transform(source_system, target_system, use_alpha):
    """Decorator for transforms."""
    def decorator(transform_function):
        transform_parameters = signature(transform_function).parameters
        uses_out_angles = "out_angles" in transform_parameters

        def wrapper(input_collection, *args, **kwargs):
            if uses_out_angles and "out_angles" not in transform_parameters:
                msg = "Out angles is expected but not provided for this function"
                raise InvalidDataError(msg)
            if use_alpha and "alpha" not in input_collection:
                msg = "alpha expected in input_collection but not found."
                raise MissingAlphaError(msg)
            required_keys = SYSTEM_REQUIRED_KEYS[source_system]
            for key in required_keys:
                if key not in input_collection:
                    msg = f"Expected key of {key} for {source_system} but not found."
                    raise InvalidDataError(msg)
            return transform_function(input_collection, *args, **kwargs)

        wrapper.uses_out_angles = uses_out_angles
        wrapper.uses_alpha = use_alpha
        return wrapper
    return decorator


@transform(System.npol, System.mzpsolar, use_alpha=False)
@u.quantity_input
def npol_to_mzpsolar(input_collection, reference_angle=0 * u.degree, **kwargs):
    """
    Notes
    ------
    Equation 44 in DeForest et al. 2022.
    """
    input_keys = list(input_collection.keys())
    phi = [input_collection[key].meta['POLAR'] for key in input_keys if key != 'alpha'] * u.degree
    mzp_angles = [-60, 0, 60] * u.degree  # theta angle in Eq 44

    data_shape = input_collection[input_keys[0]].data.shape
    data_npol = np.zeros([data_shape[0], data_shape[1], 3, 1])

    conv_matrix = np.array([[(4 * np.cos(phi[i] - mzp_angles[j] - reference_angle) ** 2 - 1) / 3
                             for j in range(3)] for i in range(3)])

    for i, key in enumerate(key for key in input_keys if key != 'alpha'):
        data_npol[:, :, i, 0] = input_collection[key].data

    try:
        conv_matrix_inv = np.linalg.inv(conv_matrix)
    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            raise SolpolpyError("Conversion matrix is degenerate")

    data_mzp_solar = np.matmul(conv_matrix_inv, data_npol)

    metas = [copy.copy(input_collection[input_keys[0]].meta) for _ in range(3)]

    for meta, original_angle, target_angle in zip(metas, input_keys, [-60, 0, 60] * u.degree):
        meta.update({
            'POLAR': target_angle,
            'POLARREF': "Solar",
            "POLAROFF": input_collection[original_angle].meta.get("POLAROFF", 0 * u.degree)
        })

    mask = combine_all_collection_masks(input_collection)
    cube_list = [(key, NDCube(data_mzp_solar[:, :, i, 0], wcs=input_collection[input_keys[0]].wcs,
                              mask=mask, meta=metas[i])) for i, key in enumerate(["M", "Z", "P"])]

    for p_angle in input_keys:
        if p_angle.lower() == "alpha":
            cube_list.append(("alpha", NDCube(input_collection["alpha"].data * u.radian,
                                              wcs=input_collection[input_keys[0]].wcs,
                                              meta=input_collection["alpha"].meta)))
    return NDCollection(cube_list, meta={}, aligned_axes="all")


@transform(System.mzpsolar, System.bpb, use_alpha=True)
def mzpsolar_to_bpb(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 7 and 9 in DeForest et al. 2022.

    """""
    # TODO: need to check if 3 angles are input.
    # TODO: need to check if separated appropriately if not create quality warning.
    input_dict = {}
    in_list = list(input_collection)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[input_collection[p_angle].meta["POLAR"]] = input_collection[p_angle].data

    alpha = input_collection["alpha"].data * u.radian
    B = (2 / 3) * (np.sum([ith_polarizer_brightness
                           for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))
    metaB, metapB = copy.copy(input_collection["M"].meta), copy.copy(input_collection["M"].meta)
    metaB.update(POLAR="B")
    metapB.update(POLAR="pB")

    mask = combine_all_collection_masks(input_collection)

    BpB_cube = [("B", NDCube(B, wcs=input_collection["M"].wcs, mask=mask, meta=metaB)),
                ("pB", NDCube(pB, wcs=input_collection["M"].wcs, mask=mask, meta=metapB)),
                ("alpha", NDCube(alpha, wcs=input_collection["M"].wcs, mask=mask))]
    # TODO: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


@transform(System.bpb, System.mzpsolar, use_alpha=True)
def bpb_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 4 in DeForest et al. 2022.
    """
    alpha = input_collection["alpha"].data * u.radian
    B, pB = input_collection["B"].data, input_collection["pB"].data
    mzp_angles = [-60, 0, 60] * u.degree
    Bmzp = {}
    for angle in mzp_angles:
        Bmzp[angle] = (1 / 2) * (B - pB * (np.cos(2 * (angle - alpha))))

    metaM, metaZ, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["B"].meta), copy.copy(
        input_collection["B"].meta)
    metaM.update(POLAR=-60*u.degree, POLARREF='Solar')
    metaZ.update(POLAR=0*u.degree, POLARREF='Solar')
    metaP.update(POLAR=60*u.degree, POLARREF='Solar')
    mask = combine_all_collection_masks(input_collection)
    Bmzp_cube = [("M", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaM)),
                 ("Z", NDCube(Bmzp[0 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaZ)),
                 ("P", NDCube(Bmzp[60 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaP)),
                 ("alpha", NDCube(alpha, wcs=input_collection["B"].wcs))]
    # TODO: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


@transform(System.bpb, System.btbr, use_alpha=True)
def bpb_to_btbr(input_collection, **kwargs):
    """Notes
    -----
    Equation 1 and 2 in DeForest et al. 2022.
    """
    alpha = input_collection["alpha"].data * u.radian
    B, pB = input_collection["B"].data, input_collection["pB"].data
    Br = (B - pB) / 2
    Bt = (B + pB) / 2

    metaBr, metaBt = copy.copy(input_collection["B"].meta), copy.copy(input_collection["B"].meta)
    metaBr.update(POLAR="Br")
    metaBt.update(POLAR="Bt")

    mask = combine_all_collection_masks(input_collection)
    BtBr_cube = [("Bt", NDCube(Bt, wcs=input_collection["B"].wcs, mask=mask, meta=metaBt)),
                 ("Br", NDCube(Br, wcs=input_collection["B"].wcs, mask=mask, meta=metaBr)),
                 ("alpha", NDCube(alpha, wcs=input_collection["B"].wcs))]
    # TODO: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BtBr_cube, meta={}, aligned_axes="all")


@transform(System.btbr, System.bpb, use_alpha=True)
def btbr_to_bpb(input_collection, **kwargs):
    """Notes
    -----
    Equation in Table 1 in DeForest et al. 2022.
    """
    if "alpha" not in input_collection:
        msg = "missing alpha"
        raise ValueError(msg)

    alpha = input_collection["alpha"].data * u.radian
    Bt, Br = input_collection["Bt"].data, input_collection["Br"].data
    pB = (Bt - Br)
    B = (Bt + Br)

    metaB, metapB = copy.copy(input_collection["Bt"].meta), copy.copy(input_collection["Bt"].meta)
    metaB.update(POLAR="B")
    metapB.update(POLAR="pB")

    mask = combine_all_collection_masks(input_collection)
    BpB_cube = [("B", NDCube(B, wcs=input_collection["Bt"].wcs, mask=mask, meta=metaB)),
                ("pB", NDCube(pB, wcs=input_collection["Bt"].wcs, mask=mask, meta=metapB)),
                ("alpha", NDCube(alpha, wcs=input_collection["Bt"].wcs, mask=mask))]
    # TODO: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(BpB_cube, meta={}, aligned_axes="all")


@transform(System.mzpsolar, System.stokes, use_alpha=False)
def mzpsolar_to_stokes(input_collection, **kwargs):
    """Notes
    -----
    Equation 9, 12 and 13 in DeForest et al. 2022.
    """
    Bm, Bz, Bp = input_collection["M"].data, input_collection["Z"].data, input_collection["P"].data

    mueller_matrix = (2 / 3) * np.array([[1, 1, 1], [-1, 2, -1], [-np.sqrt(3), 0, np.sqrt(3)]])

    Bi = mueller_matrix[0, 0] * Bm + mueller_matrix[0, 1] * Bz + mueller_matrix[0, 2] * Bp
    Bq = mueller_matrix[1, 0] * Bm + mueller_matrix[1, 1] * Bz + mueller_matrix[1, 2] * Bp
    Bu = mueller_matrix[2, 0] * Bm + mueller_matrix[2, 1] * Bz + mueller_matrix[2, 2] * Bp

    metaI, metaQ, metaU = copy.copy(input_collection["M"].meta), copy.copy(input_collection["Z"].meta), copy.copy(
        input_collection["P"].meta)
    metaI.update(POLAR="Stokes I")
    metaQ.update(POLAR="Stokes Q")
    metaU.update(POLAR="Stokes U")

    mask = combine_all_collection_masks(input_collection)
    BStokes_cube = [("I", NDCube(Bi, wcs=input_collection["M"].wcs, mask=mask, meta=metaI)),
                    ("Q", NDCube(Bq, wcs=input_collection["M"].wcs, mask=mask, meta=metaQ)),
                    ("U", NDCube(Bu, wcs=input_collection["M"].wcs, mask=mask, meta=metaU))]
    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")


@transform(System.stokes, System.mzpsolar, use_alpha=False)
def stokes_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 11 in DeForest et al. 2022. with alpha = np.pi/2
    """
    alpha = np.pi / 2
    Bi, Bq, Bu = input_collection["I"].data, input_collection["Q"].data, input_collection["U"].data

    inv_mul_mx = (1 / 2) * np.array([[1, -np.cos(2 * (-np.pi / 3 - alpha)), -np.sin(2 * (-np.pi / 3 - alpha))],
                                     [1, -np.cos(2 * (0 - alpha)), 0],
                                     [1, -np.cos(2 * (np.pi / 3 - alpha)), -np.sin(2 * (np.pi / 3 - alpha))]])

    Bm = inv_mul_mx[0, 0] * Bi + inv_mul_mx[0, 1] * Bq + inv_mul_mx[0, 2] * Bu
    Bz = inv_mul_mx[1, 0] * Bi + inv_mul_mx[1, 1] * Bq + inv_mul_mx[1, 2] * Bu
    Bp = inv_mul_mx[2, 0] * Bi + inv_mul_mx[2, 1] * Bq + inv_mul_mx[2, 2] * Bu

    metaM, metaZ, metaP = copy.copy(input_collection["I"].meta), copy.copy(input_collection["Q"].meta), copy.copy(
        input_collection["U"].meta)
    metaM.update(POLAR=-60*u.degree, POLARREF='Solar')
    metaZ.update(POLAR=0*u.degree, POLARREF='Solar')
    metaP.update(POLAR=60*u.degree, POLARREF='Solar')

    mask = combine_all_collection_masks(input_collection)
    Bmzp_cube = [("M", NDCube(Bm, wcs=input_collection["I"].wcs, mask=mask, meta=metaM)),
                 ("Z", NDCube(Bz, wcs=input_collection["I"].wcs, mask=mask, meta=metaZ)),
                 ("P", NDCube(Bp, wcs=input_collection["I"].wcs, mask=mask, meta=metaP)),
                 ("alpha", NDCube(np.zeros_like(Bm) + alpha, wcs=input_collection["I"].wcs, mask=mask))]

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


@transform(System.mzpsolar, System.bp3, use_alpha=True)
def mzpsolar_to_bp3(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 7, 9 and 10 in DeForest et al. 2022.
    """""
    input_dict = {}
    in_list = list(input_collection)

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[input_collection[p_angle].meta["POLAR"]] = input_collection[p_angle].data

    alpha = input_collection["alpha"].data * u.radian
    B = (2 / 3) * (np.sum([ith_polarizer_brightness for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))

    pBp = (-4 / 3) * (np.sum([ith_polarizer_brightness * np.sin(2 * (ith_angle - alpha))
                              for ith_angle, ith_polarizer_brightness
                              in input_dict.items() if ith_angle != "alpha"], axis=0))
    # TODO: update header properly
    metaB, metapB, metapBp = copy.copy(input_collection["M"].meta), copy.copy(input_collection["M"].meta), copy.copy(
        input_collection["M"].meta)
    metaB.update(POLAR="B")
    metapB.update(POLAR="pB")
    metapBp.update(POLAR="pB-prime")

    mask = combine_all_collection_masks(input_collection)
    Bp3_cube = [("B", NDCube(B, wcs=input_collection["M"].wcs, mask=mask, meta=metaB)),
                ("pB", NDCube(pB, wcs=input_collection["M"].wcs, mask=mask, meta=metapB)),
                ("pBp", NDCube(pBp, wcs=input_collection["M"].wcs, mask=mask, meta=metapBp)),
                ("alpha", NDCube(alpha, wcs=input_collection["M"].wcs, mask=mask))]
    # TODO: WCS for alpha needs to be generated wrt to solar north

    return NDCollection(Bp3_cube, meta={}, aligned_axes="all")


@transform(System.bp3, System.mzpsolar, use_alpha=True)
def bp3_to_mzpsolar(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022.
    """""
    B, pB, pBp = input_collection["B"].data, input_collection["pB"].data, input_collection["pBp"].data
    alpha = input_collection["alpha"].data * u.radian

    mzp_angles = [-60, 0, 60] * u.degree
    Bmzp = {}
    for angle in mzp_angles:
        Bmzp[angle] = (1 / 2) * (B - np.cos(2 * (angle - alpha)) * pB -
                               np.cos(2 * (angle - alpha)) * pBp)

    metaM, metaZ, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta), copy.copy(
        input_collection["pBp"].meta)
    metaM.update(POLAR=-60*u.degree, POLARREF='Solar')
    metaZ.update(POLAR=0*u.degree, POLARREF='Solar')
    metaP.update(POLAR=60*u.degree, POLARREF='Solar')

    mask = combine_all_collection_masks(input_collection)
    Bmzp_cube = [("M", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaM)),
                 ("Z", NDCube(Bmzp[0 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaZ)),
                 ("P", NDCube(Bmzp[60 * u.degree], wcs=input_collection["B"].wcs, mask=mask, meta=metaP)),
                 ("alpha", NDCube(alpha, wcs=input_collection["B"].wcs, mask=mask))]

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


@transform(System.btbr, System.mzpsolar, use_alpha=True)
def btbr_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 3 in DeForest et al. 2022.
    """
    alpha = input_collection["alpha"].data * u.radian
    Bt = input_collection["Bt"].data
    Br = input_collection["Br"].data

    mzp_angles = [-60, 0, 60] * u.degree
    Bmzp = {}
    for angle in mzp_angles:
        Bmzp[angle] = Bt * (np.sin(angle - alpha)) ** 2 + Br * (np.cos(angle - alpha)) ** 2

    metaM, metaZ, metaP = copy.copy(input_collection["Bt"].meta), copy.copy(input_collection["Bt"].meta), copy.copy(
        input_collection["Bt"].meta)
    metaM.update(POLAR=-60*u.degree, POLARREF='Solar')
    metaZ.update(POLAR=0*u.degree, POLARREF='Solar')
    metaP.update(POLAR=60*u.degree, POLARREF='Solar')

    mask = combine_all_collection_masks(input_collection)
    Bmzp_cube = [("M", NDCube(Bmzp[-60 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask, meta=metaM)),
                 ("Z", NDCube(Bmzp[0 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask, meta=metaZ)),
                 ("P", NDCube(Bmzp[60 * u.degree], wcs=input_collection["Bt"].wcs, mask=mask, meta=metaP)),
                 ("alpha", NDCube(alpha, wcs=input_collection["Bt"].wcs, mask=mask))]

    return NDCollection(Bmzp_cube, meta={}, aligned_axes="all")


@transform(System.bp3, System.bthp, use_alpha=True)
def bp3_to_bthp(input_collection, **kwargs):
    """
    Notes
    ------
    Equations 9, 15, 16 in DeForest et al. 2022.
    """""
    B, pB, pBp = input_collection["B"].data, input_collection["pB"].data, input_collection["pBp"].data
    alpha = input_collection["alpha"].data * u.radian

    theta = (1 / 2) * np.arctan2(pBp, pB) * u.radian + np.pi / 2 * u.radian + alpha
    p = np.sqrt(pB ** 2 + pBp ** 2) / B

    metaTh, metaP = copy.copy(input_collection["B"].meta), copy.copy(input_collection["pB"].meta)
    metaTh.update(POLAR="Theta")
    metaP.update(POLAR="Degree of Polarization")

    mask = combine_all_collection_masks(input_collection)
    Bthp_cube = [("B", NDCube(B, wcs=input_collection["B"].wcs, mask=mask, meta=input_collection["B"].meta)),
                 ("theta", NDCube(theta, wcs=input_collection["B"].wcs, mask=mask, meta=metaTh)),
                 ("p", NDCube(p, wcs=input_collection["B"].wcs, mask=mask, meta=metaP))]

    return NDCollection(Bthp_cube, meta={}, aligned_axes="all")


@transform(System.btbr, System.npol, use_alpha=True)
@u.quantity_input
def btbr_to_npol(input_collection, out_angles: u.degree, **kwargs):
    """Notes
    -----
    Equation 3 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    alpha = input_collection["alpha"].data * u.radian
    Bt, Br = input_collection["Bt"].data, input_collection["Br"].data

    Bnpol = {}
    Bnpol_cube = []
    mask = combine_all_collection_masks(input_collection)
    for angle in out_angles:
        Bnpol[angle] = Bt * (np.sin(angle - alpha)) ** 2 + Br * (np.cos(angle - alpha)) ** 2
        meta_tmp = copy.copy(input_collection["Bt"].meta)
        meta_tmp.update(POLAR=angle)
        Bnpol_cube.append((str(angle), NDCube(Bnpol[angle], wcs=input_collection["Bt"].wcs, mask=mask,  meta=meta_tmp)))
    Bnpol_cube.append(("alpha", NDCube(alpha, wcs=input_collection["Bt"].wcs, mask=mask)))

    return NDCollection(Bnpol_cube, meta={}, aligned_axes="all")


@transform(System.mzpsolar, System.npol, use_alpha=False)
@u.quantity_input
def mzpsolar_to_npol(input_collection, out_angles: u.degree, reference_angle=0 * u.degree, **kwargs):
    """Notes
    -----
    Equation 45 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    in_list = list(input_collection)
    input_dict = {}

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[input_collection[p_angle].meta["POLAR"]] = input_collection[p_angle].data

    output_cubes = []
    mask = combine_all_collection_masks(input_collection)
    first_meta = input_collection[in_list[0]].meta
    first_wcs = input_collection[in_list[0]].wcs
    for out_angle in out_angles:
        value = (1/3) * np.sum([input_cube.data * (4 * np.square(np.cos(out_angle - input_angle - reference_angle)) - 1)
                                     for input_angle, input_cube in input_dict.items()], axis=0)
        out_meta = copy.copy(first_meta)
        out_meta.update(POLAR=out_angle)
        output_cubes.append((str(out_angle),
                           NDCube(value, wcs=first_wcs, mask=mask, meta=out_meta)))

    if "alpha" in input_collection:
        alpha = input_collection["alpha"].data * u.radian
        output_cubes.append(("alpha", NDCube(alpha, wcs=input_collection[in_list[0]].wcs, mask=mask)))

    return NDCollection(output_cubes, meta={}, aligned_axes="all")


@transform(System.fourpol, System.stokes, use_alpha=False)
def fourpol_to_stokes(input_collection, **kwargs):
    """
    Notes
    ------
    Table 1 in DeForest et al. 2022.

    """""
    Bi = input_collection[str(0 * u.degree)].data + input_collection[str(90 * u.degree)].data
    Bq = input_collection[str(90 * u.degree)].data - input_collection[str(0 * u.degree)].data
    Bu = input_collection[str(135 * u.degree)].data - input_collection[str(45 * u.degree)].data

    metaI, metaQ, metaU = (copy.copy(input_collection[str(0 * u.degree)].meta),
                           copy.copy(input_collection[str(0 * u.degree)].meta),
                           copy.copy(input_collection[str(0 * u.degree)].meta))
    metaI.update(POLAR="Stokes I")
    metaQ.update(POLAR="Stokes Q")
    metaU.update(POLAR="Stokes U")

    mask = combine_all_collection_masks(input_collection)
    BStokes_cube = [("I", NDCube(Bi, wcs=input_collection[str(0 * u.degree)].wcs, mask=mask, meta=metaI)),
                    ("Q", NDCube(Bq, wcs=input_collection[str(0 * u.degree)].wcs, mask=mask, meta=metaQ)),
                    ("U", NDCube(Bu, wcs=input_collection[str(0 * u.degree)].wcs, mask=mask, meta=metaU))]

    return NDCollection(BStokes_cube, meta={}, aligned_axes="all")


@transform(System.mzpsolar, System.mzpinstru, use_alpha=False)
def mzpsolar_to_mzpinstru(input_collection, reference_angle=0 * u.degree, **kwargs):
    """Notes
        -----
        Equation 45 in DeForest et al. 2022.
        out_angles: list of target angles in degree
        """
    in_list = list(input_collection)
    input_dict = {}

    for p_angle in in_list:
        if p_angle == "alpha":
            break
        input_dict[input_collection[p_angle].meta["POLAR"]] = input_collection[p_angle].data

    output_cubes = []
    mask = combine_all_collection_masks(input_collection)
    satellite_orientation = extract_crota_from_wcs(input_collection['Z'])
    mzp_angles = [input_collection[key].meta['POLAR'] for key in list(input_collection.keys()) if key != 'alpha']*u.degree
    out_angles = mzp_angles + satellite_orientation

    for out_angle, key in zip(out_angles, ["M", "Z", "P"]):
        value = (1 / 3) * np.sum(
            [input_cube.data * (4 * np.square(np.cos(out_angle - input_angle - reference_angle)) - 1)
             for input_angle, input_cube in input_dict.items()], axis=0)
        out_meta = copy.copy(input_collection[key].meta)
        out_meta.update(POLARREF="Instrument")
        output_cubes.append((key,
                             NDCube(value, wcs=input_collection[key].wcs, mask=mask, meta=out_meta)))

    if "alpha" in input_collection:
        alpha = input_collection["alpha"].data * u.radian
        output_cubes.append(("alpha", NDCube(alpha, wcs=input_collection[in_list[0]].wcs, mask=mask)))

    return NDCollection(output_cubes, meta={}, aligned_axes="all")


@transform(System.mzpinstru, System.mzpsolar, use_alpha=False)
def mzpinstru_to_mzpsolar(input_collection, reference_angle=0*u.degree, **kwargs):
    """
    Notes
    ------
    For rotating frames like NFI, ASPIICS, CODEX
    Input has MZP in instrument reference.
    Equation 44 in DeForest et al. 2022.
    """
    input_keys = list(input_collection.keys())
    satellite_orientation = extract_crota_from_wcs(input_collection['Z'])
    polarizer_difference = [input_collection[k].meta['POLAROFF'] if 'POLAROFF' in input_collection[k].meta else 0
                            for k in ['M', 'Z', 'P']] * u.degree
    phi = [input_collection[key].meta['POLAR'] for key in input_keys if key != 'alpha']*u.degree
    phi = phi + satellite_orientation + polarizer_difference
    mzp_angles = [-60, 0, 60] * u.degree    # theta angle in Eq 44
    data_shape = input_collection[input_keys[0]].data.shape
    data_npol = np.zeros([data_shape[0], data_shape[1], 3, 1])
    conv_matrix = np.array([[(4 * np.cos(phi[i] - mzp_angles[j] - reference_angle) ** 2 - 1) / 3
                             for j in range(3)] for i in range(3)])

    for i, key in enumerate(key for key in input_keys if key != 'alpha'):
        data_npol[:, :, i, 0] = input_collection[key].data
    try:
        conv_matrix_inv = np.linalg.inv(conv_matrix)
    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            raise SolpolpyError("Conversion matrix is degenerate")

    data_mzp_solar = np.matmul(conv_matrix_inv, data_npol)
    mask = combine_all_collection_masks(input_collection)

    metas = [{'POLAR': target_angle,
              'POLARREF': "Solar",
              "POLAROFF": input_collection[original_angle].meta.get("POLAROFF", 0*u.degree)}
             for original_angle, target_angle in zip(input_keys, [-60, 0, 60] * u.degree)]

    cube_list = [(key, NDCube(data_mzp_solar[:, :, i, 0], wcs=input_collection[input_keys[0]].wcs,
                mask=mask, meta=metas[i])) for i, key in enumerate(["M", "Z", "P"])]
    for p_angle in input_keys:
        if p_angle.lower() == "alpha":
            cube_list.append(("alpha", NDCube(input_collection["alpha"].data * u.radian,
                                              wcs=input_collection[input_keys[0]].wcs,
                                              meta=input_collection["alpha"].meta)))
    return NDCollection(cube_list, meta={}, aligned_axes="all")

# Build the graph at the bottom so all transforms are defined
transform_graph = nx.DiGraph()
transform_functions = getmembers(sys.modules[__name__], isfunction)
for function_name, function in transform_functions:
    if "_to_" in function_name:
        source, destination = function_name.split("_to_")
        transform_graph.add_edge(System[source.lower()],
                                 System[destination.lower()],
                                 func=function)
