"""Polarizer functions for solpolpy."""

import sys
import warnings
from enum import StrEnum
from inspect import signature, getmembers, isfunction

import astropy.units as u
import networkx as nx
import numpy as np
from ndcube import NDCollection, NDCube

from solpolpy.errors import InvalidDataError, MissingAlphaError
from solpolpy.physics import (
    MZP_ANGLES,
    angle_difference_radians,
    as_angle,
    bp3_from_analyzer_brightness,
    bp3_to_analyzer_brightness,
    clone_meta,
    combine_mask,
    data_keys,
    project_three_polarizer_brightness,
    solve_three_polarizer_brightness,
    stack_data,
    template_cube,
    wrap_linear_polarization,
)
from solpolpy.util import compute_lats, solnorth_from_wcs

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
        transform_signature = signature(transform_function)
        transform_parameters = transform_signature.parameters
        uses_out_angles = "out_angles" in transform_parameters
        uses_in_angles = "in_angles" in transform_parameters
        requires_out_angles = uses_out_angles and transform_parameters["out_angles"].default is transform_signature.empty
        requires_in_angles = uses_in_angles and transform_parameters["in_angles"].default is transform_signature.empty

        def wrapper(input_collection, *args, **kwargs):
            bound_args = transform_signature.bind_partial(input_collection, *args, **kwargs)
            if requires_out_angles and bound_args.arguments.get("out_angles") is None:
                msg = "Out angles is expected but not provided for this function"
                raise InvalidDataError(msg)
            if requires_in_angles and bound_args.arguments.get("in_angles") is None:
                msg = "In angles is expected but not provided for this function"
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
        wrapper.uses_in_angles = uses_in_angles
        wrapper.requires_out_angles = requires_out_angles
        wrapper.requires_in_angles = requires_in_angles
        wrapper.uses_alpha = use_alpha
        wrapper.fcn = transform_function
        return wrapper
    return decorator


def _alpha_data(collection: NDCollection):
    return as_angle(collection["alpha"].data, u.radian).to(u.radian)


def _collection_from_cubes(cubes):
    return NDCollection(cubes, meta={}, aligned_axes="all")


def _shared_polaroff(input_collection: NDCollection, keys):
    offsets = []
    for key in keys:
        if key in input_collection and "POLAROFF" in input_collection[key].meta:
            offsets.append(as_angle(input_collection[key].meta["POLAROFF"], u.degree))

    if not offsets:
        return None

    first = offsets[0].to_value(u.degree)
    if all(np.allclose(offset.to_value(u.degree), first) for offset in offsets[1:]):
        return offsets[0]

    warnings.warn(
        "Input collection contains mixed POLAROFF values; reduced output metadata cannot preserve per-channel offsets.",
        stacklevel=2,
    )
    return None


def _meta_with_shared_polaroff(cube, shared_polaroff, **updates):
    meta = clone_meta(cube, **updates)
    if shared_polaroff is None:
        meta.pop("POLAROFF", None)
    else:
        meta["POLAROFF"] = shared_polaroff
    return meta


def _warn_if_information_is_lost(pBp, B, context):
    with np.errstate(divide="ignore", invalid="ignore"):
        relative_pbp = np.abs(
            np.divide(
                pBp,
                B,
                out=np.zeros_like(np.asarray(pBp, dtype=float), dtype=float),
                where=np.not_equal(B, 0),
            )
        )

    max_relative_pbp = np.nanmax(relative_pbp)
    if max_relative_pbp > 1e-6:
        warnings.warn(
            (
                f"{context} assumes pBp = 0, so roundtrip recovery will not be exact when the input carries "
                f"nonzero pBp. max(|pBp/B|)={max_relative_pbp:.3%}"
            ),
            stacklevel=2,
        )


def _append_alpha(cubes, input_collection: NDCollection, mask=None):
    if "alpha" not in input_collection:
        return cubes
    alpha_cube = input_collection["alpha"]
    cubes.append(
        (
            "alpha",
            NDCube(
                _alpha_data(input_collection),
                wcs=alpha_cube.wcs,
                mask=mask if mask is not None else alpha_cube.mask,
                meta=alpha_cube.meta,
            ),
        )
    )
    return cubes


def _mzp_cubes_from_stack(data_stack, input_collection: NDCollection, mask=None, preferred_key: str | None = None):
    cube_template = template_cube(input_collection, preferred_key=preferred_key)
    shared_polaroff = _shared_polaroff(input_collection, data_keys(input_collection))
    cubes = []
    for key, angle, data in zip(["M", "Z", "P"], MZP_ANGLES, data_stack, strict=False):
        cubes.append(
            (
                key,
                NDCube(
                    data,
                    wcs=cube_template.wcs,
                    mask=mask,
                    meta=_meta_with_shared_polaroff(cube_template, shared_polaroff, POLAR=angle, POLARREF="Solar"),
                ),
            )
        )
    return cubes


def _instrument_frame_analyzer_angles(input_collection: NDCollection):
    data_shape = template_cube(input_collection, preferred_key="Z").data.shape
    lats = compute_lats(input_collection["Z"].wcs, data_shape)

    polarizer_difference = {
        key: as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree) for key in ["M", "Z", "P"]
    }
    solar_north = {
        key: solnorth_from_wcs(input_collection[key].wcs, shape=data_shape, precomputed_lats=lats)
        for key in ["M", "Z", "P"]
    }

    return np.stack(
        [
            wrap_linear_polarization(solar_north["M"] + MZP_ANGLES[0] + polarizer_difference["M"]),
            wrap_linear_polarization(solar_north["Z"] + MZP_ANGLES[1] + polarizer_difference["Z"]),
            wrap_linear_polarization(solar_north["P"] + MZP_ANGLES[2] + polarizer_difference["P"]),
        ]
    )


def _solar_north_angles(input_collection: NDCollection):
    data_shape = template_cube(input_collection, preferred_key="Z").data.shape
    lats = compute_lats(input_collection["Z"].wcs, data_shape)
    return {
        key: solnorth_from_wcs(input_collection[key].wcs, shape=data_shape, precomputed_lats=lats)
        for key in ["M", "Z", "P"]
    }


@transform(System.mzpsolar, System.bpb, use_alpha=True)
def mzpsolar_to_bpb(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 7 and 9 in DeForest et al. 2022.

    """""
    alpha = _alpha_data(input_collection)
    analyzer_stack = np.stack(stack_data(input_collection, ["M", "Z", "P"]), axis=0)
    B, pB, pBp = bp3_from_analyzer_brightness(analyzer_stack, MZP_ANGLES, alpha)
    _warn_if_information_is_lost(pBp, B, "mzpsolar_to_bpb")

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="M")
    shared_polaroff = _shared_polaroff(input_collection, ["M", "Z", "P"])
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="B"))),
        ("pB", NDCube(pB, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="pB"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bpb, System.mzpsolar, use_alpha=True)
def bpb_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 4 in DeForest et al. 2022.
    """
    alpha = _alpha_data(input_collection)
    B, pB = input_collection["B"].data, input_collection["pB"].data
    mzp_stack = bp3_to_analyzer_brightness(B, pB, np.zeros_like(pB), alpha, MZP_ANGLES)
    mask = combine_mask(input_collection)
    cubes = _mzp_cubes_from_stack(mzp_stack, input_collection, mask=mask, preferred_key="B")
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bpb, System.btbr, use_alpha=True)
def bpb_to_btbr(input_collection, **kwargs):
    """Notes
    -----
    Equation 1 and 2 in DeForest et al. 2022.
    """
    alpha = _alpha_data(input_collection)
    B, pB = input_collection["B"].data, input_collection["pB"].data
    Br = (B - pB) / 2
    Bt = (B + pB) / 2

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="B")
    cubes = [
        ("Bt", NDCube(Bt, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Bt"))),
        ("Br", NDCube(Br, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Br"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.btbr, System.bpb, use_alpha=True)
def btbr_to_bpb(input_collection, **kwargs):
    """Notes
    -----
    Equation in Table 1 in DeForest et al. 2022.
    """
    alpha = _alpha_data(input_collection)
    Bt, Br = input_collection["Bt"].data, input_collection["Br"].data
    pB = (Bt - Br)
    B = (Bt + Br)

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="Bt")
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="B"))),
        ("pB", NDCube(pB, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="pB"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


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

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="M")
    shared_polaroff = _shared_polaroff(input_collection, ["M", "Z", "P"])
    cubes = [
        ("I", NDCube(Bi, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="Stokes I"))),
        ("Q", NDCube(Bq, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="Stokes Q"))),
        ("U", NDCube(Bu, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="Stokes U"))),
    ]
    return _collection_from_cubes(cubes)


@transform(System.stokes, System.mzpsolar, use_alpha=False)
def stokes_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 11 in DeForest et al. 2022. with alpha = np.pi/2
    """
    alpha = 90 * u.degree
    Bi, Bq, Bu = input_collection["I"].data, input_collection["Q"].data, input_collection["U"].data

    inv_mul_mx = (1 / 2) * np.array([
        [1, -np.cos(2 * angle_difference_radians(-60 * u.degree, alpha)), -np.sin(2 * angle_difference_radians(-60 * u.degree, alpha))],
        [1, -np.cos(2 * angle_difference_radians(0 * u.degree, alpha)), 0],
        [1, -np.cos(2 * angle_difference_radians(60 * u.degree, alpha)), -np.sin(2 * angle_difference_radians(60 * u.degree, alpha))],
    ])

    Bm = inv_mul_mx[0, 0] * Bi + inv_mul_mx[0, 1] * Bq + inv_mul_mx[0, 2] * Bu
    Bz = inv_mul_mx[1, 0] * Bi + inv_mul_mx[1, 1] * Bq + inv_mul_mx[1, 2] * Bu
    Bp = inv_mul_mx[2, 0] * Bi + inv_mul_mx[2, 1] * Bq + inv_mul_mx[2, 2] * Bu

    mask = combine_mask(input_collection)
    cubes = _mzp_cubes_from_stack([Bm, Bz, Bp], input_collection, mask=mask, preferred_key="I")
    alpha_plane = np.full(np.shape(Bm), alpha.to_value(u.radian)) * u.radian
    cubes.append(("alpha", NDCube(alpha_plane, wcs=input_collection["I"].wcs, mask=mask)))
    return _collection_from_cubes(cubes)


@transform(System.mzpsolar, System.bp3, use_alpha=True)
def mzpsolar_to_bp3(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 7, 9 and 10 in DeForest et al. 2022.
    """""
    alpha = _alpha_data(input_collection)
    analyzer_stack = np.stack(stack_data(input_collection, ["M", "Z", "P"]), axis=0)
    B, pB, pBp = bp3_from_analyzer_brightness(analyzer_stack, MZP_ANGLES, alpha)

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="M")
    shared_polaroff = _shared_polaroff(input_collection, ["M", "Z", "P"])
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="B"))),
        ("pB", NDCube(pB, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="pB"))),
        ("pBp", NDCube(pBp, wcs=template.wcs, mask=mask, meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR="pB-prime"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bp3, System.mzpsolar, use_alpha=True)
def bp3_to_mzpsolar(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 11 in DeForest et al. 2022.
    """""
    B, pB, pBp = input_collection["B"].data, input_collection["pB"].data, input_collection["pBp"].data
    alpha = _alpha_data(input_collection)
    mzp_stack = bp3_to_analyzer_brightness(B, pB, pBp, alpha, MZP_ANGLES)

    mask = combine_mask(input_collection)
    cubes = _mzp_cubes_from_stack(mzp_stack, input_collection, mask=mask, preferred_key="B")
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.btbr, System.mzpsolar, use_alpha=True)
def btbr_to_mzpsolar(input_collection, **kwargs):
    """Notes
    -----
    Equation 3 in DeForest et al. 2022.
    """
    alpha = _alpha_data(input_collection)
    Bt = input_collection["Bt"].data
    Br = input_collection["Br"].data

    mzp_stack = np.stack(
        [
            Bt * np.sin(angle_difference_radians(angle, alpha)) ** 2 + Br * np.cos(angle_difference_radians(angle, alpha)) ** 2
            for angle in MZP_ANGLES
        ],
        axis=0,
    )
    mask = combine_mask(input_collection)
    cubes = _mzp_cubes_from_stack(mzp_stack, input_collection, mask=mask, preferred_key="Bt")
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bp3, System.bthp, use_alpha=True)
def bp3_to_bthp(input_collection, **kwargs):
    """
    Notes
    ------
    Equations 9, 15, 16 in DeForest et al. 2022.
    """""
    B, pB, pBp = input_collection["B"].data, input_collection["pB"].data, input_collection["pBp"].data
    alpha = _alpha_data(input_collection)

    theta = wrap_linear_polarization(0.5 * np.arctan2(pBp, pB) * u.radian + np.pi / 2 * u.radian + alpha)
    p = np.divide(
        np.sqrt(pB ** 2 + pBp ** 2),
        B,
        out=np.full_like(B, np.nan, dtype=float),
        where=np.not_equal(B, 0),
    )

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="B")
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=template.meta)),
        ("theta", NDCube(theta, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Theta"))),
        ("p", NDCube(p, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Degree of Polarization"))),
    ]
    return _collection_from_cubes(cubes)


@transform(System.btbr, System.npol, use_alpha=True)
@u.quantity_input
def btbr_to_npol(input_collection, out_angles: u.degree, **kwargs):
    """Notes
    -----
    Equation 3 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    alpha = _alpha_data(input_collection)
    Bt, Br = input_collection["Bt"].data, input_collection["Br"].data

    cubes = []
    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key="Bt")
    for angle in out_angles:
        value = Bt * np.sin(angle_difference_radians(angle, alpha)) ** 2 + Br * np.cos(angle_difference_radians(angle, alpha)) ** 2
        cubes.append((str(angle), NDCube(value, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR=angle))))
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.npol, System.mzpsolar, use_alpha=False)
@u.quantity_input
def npol_to_mzpsolar(input_collection, in_angles: u.degree = None, reference_angle=0 * u.degree, **kwargs):
    """
    Notes
    ------
    Equation 44 in DeForest et al. 2022.
    """
    input_keys = data_keys(input_collection)
    phi = (
        in_angles
        if in_angles is not None
        else u.Quantity(
            [
                as_angle(input_collection[key].meta["POLAR"], u.degree)
                for key in input_keys
            ]
        )
    )

    solved_stack = solve_three_polarizer_brightness(
        np.stack(stack_data(input_collection, input_keys), axis=0),
        observed_angles=phi,
        solved_angles=MZP_ANGLES,
        reference_angle=reference_angle,
    )
    mask = combine_mask(input_collection)
    cube_list = _mzp_cubes_from_stack(solved_stack, input_collection, mask=mask)
    _append_alpha(cube_list, input_collection, mask=mask)
    return _collection_from_cubes(cube_list)


@transform(System.mzpsolar, System.npol, use_alpha=False)
@u.quantity_input
def mzpsolar_to_npol(input_collection, out_angles: u.degree, reference_angle=0 * u.degree, **kwargs):
    """Notes
    -----
    Equation 45 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    in_keys = data_keys(input_collection)
    source_angles = u.Quantity(
        [
            as_angle(input_collection[key].meta["POLAR"], u.degree)
            for key in in_keys
        ]
    )
    source_stack = np.stack(stack_data(input_collection, in_keys), axis=0)
    projected = project_three_polarizer_brightness(
        source_stack,
        source_angles=source_angles,
        target_angles=out_angles,
        reference_angle=reference_angle,
    )

    mask = combine_mask(input_collection)
    template = template_cube(input_collection)
    shared_polaroff = _shared_polaroff(input_collection, in_keys)
    output_cubes = []
    for angle, value in zip(out_angles, projected, strict=False):
        out_key = str(np.round(np.mean(angle).value)) if getattr(out_angles, "ndim", 1) > 1 else str(angle)
        output_cubes.append(
            (
                out_key,
                NDCube(
                    value,
                    wcs=template.wcs,
                    mask=mask,
                    meta=_meta_with_shared_polaroff(template, shared_polaroff, POLAR=angle),
                ),
            )
        )

    _append_alpha(output_cubes, input_collection, mask=mask)
    return _collection_from_cubes(output_cubes)


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

    mask = combine_mask(input_collection)
    template = template_cube(input_collection, preferred_key=str(0 * u.degree))
    cubes = [
        ("I", NDCube(Bi, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes I"))),
        ("Q", NDCube(Bq, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes Q"))),
        ("U", NDCube(Bu, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes U"))),
    ]
    return _collection_from_cubes(cubes)


@transform(System.mzpsolar, System.mzpinstru, use_alpha=False)
def mzpsolar_to_mzpinstru(input_collection, reference_angle=0 * u.degree, **kwargs):
    """Notes
        -----
        Equation 45 in DeForest et al. 2022.
        out_angles: list of target angles in degree
        """
    mask = combine_mask(input_collection)
    solar_north = _solar_north_angles(input_collection)
    cubes = []
    input_triplet = {
        key: (as_angle(input_collection[key].meta["POLAR"], u.degree), input_collection[key].data)
        for key in ["M", "Z", "P"]
    }
    for key in ["M", "Z", "P"]:
        nominal_angle = as_angle(input_collection[key].meta["POLAR"], u.degree)
        polaroff = as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree)
        alpha_j = nominal_angle + polaroff
        phi = wrap_linear_polarization(alpha_j - solar_north[key])
        value = (1 / 3) * np.sum(
            [
                data_i * (4 * np.cos(angle_difference_radians(phi, theta_i + reference_angle)) ** 2 - 1)
                for theta_i, data_i in input_triplet.values()
            ],
            axis=0,
        )
        cubes.append(
            (
                key,
                NDCube(
                    value,
                    wcs=input_collection[key].wcs,
                    mask=mask,
                    meta=clone_meta(
                        input_collection[key],
                        POLARREF="Instrument",
                        POLAR=nominal_angle,
                        POLAROFF=polaroff,
                    ),
                ),
            )
        )
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.mzpinstru, System.mzpsolar, use_alpha=False)
def mzpinstru_to_mzpsolar(input_collection, reference_angle=0*u.degree, **kwargs):
    """Notes
    -----
    Equation 45 in DeForest et al. 2022.
    angles: list of input angles in degree
    """
    mask = combine_mask(input_collection)
    solar_north = _solar_north_angles(input_collection)
    input_triplet = {
        key: (
            as_angle(input_collection[key].meta["POLAR"], u.degree)
            + as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree),
            input_collection[key].data,
        )
        for key in ["M", "Z", "P"]
    }
    cubes = []
    for key, mzp_angle in zip(["M", "Z", "P"], MZP_ANGLES, strict=False):
        phi = wrap_linear_polarization(solar_north[key] + mzp_angle)
        value = (1 / 3) * np.sum(
            [
                data_i * (4 * np.cos(angle_difference_radians(phi, theta_i + reference_angle)) ** 2 - 1)
                for theta_i, data_i in input_triplet.values()
            ],
            axis=0,
        )
        cubes.append(
            (
                key,
                NDCube(
                    value,
                    wcs=input_collection[key].wcs,
                    mask=mask,
                    meta=clone_meta(
                        input_collection[key],
                        POLARREF="Solar",
                        POLAR=mzp_angle,
                        POLAROFF=as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree),
                    ),
                ),
            )
        )
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)

# Build the graph at the bottom so all transforms are defined
transform_graph = nx.DiGraph()
transform_functions = getmembers(sys.modules[__name__], isfunction)
for function_name, function in transform_functions:
    if "_to_" in function_name:
        source, destination = function_name.split("_to_")
        if source.lower() in System.__members__ and destination.lower() in System.__members__:
            transform_graph.add_edge(System[source.lower()],
                                     System[destination.lower()],
                                     func=function)
