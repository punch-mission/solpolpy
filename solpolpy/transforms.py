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
    as_angle,
    bp3_from_polarizer_brightness,
    bp3_to_polarizer_brightness,
    calculate_angle_difference,
    clone_meta,
    get_data_keys,
    get_template_cube,
    solve_from_three_polarizer_brightness,
    solve_to_three_polarizer_brightness,
    stack_data,
    wrap_linear_polarization,
)
from solpolpy.util import combine_all_collection_masks, compute_lats, solnorth_from_wcs

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


def transform(source_system: System, target_system: System, use_alpha: bool):
    """Decorator for transforms."""
    def decorator(transform_function):
        transform_signature = signature(transform_function)
        transform_parameters = transform_signature.parameters
        # ``uses_*`` means the transform accepts optional angle arguments that
        # may need to be forwarded through a composed path. ``requires_*`` is
        # stricter: the argument has no default and must be provided by resolve.
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


def _alpha_data(collection: NDCollection) -> u.Quantity:
    """Return the alpha cube data as radians."""
    return as_angle(collection["alpha"].data, u.radian).to(u.radian)


def _angle_difference_radians(target_angle: u.Quantity, source_angle: u.Quantity) -> np.ndarray:
    """Return a wrapped angle difference in radian values for trigonometry."""
    return calculate_angle_difference(as_angle(target_angle, u.radian), as_angle(source_angle, u.radian)).to_value(u.radian)


def _collection_from_cubes(cubes: list[tuple[str, NDCube]]) -> NDCollection:
    """Build a consistently aligned collection from transform output cubes."""
    return NDCollection(cubes, meta={}, aligned_axes="all")


def _shared_polaroff(input_collection: NDCollection, keys: list[str]) -> u.Quantity | None:
    """Return the shared POLAROFF value when output metadata can preserve it."""
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


def _meta_with_shared_polaroff(cube: NDCube, shared_polaroff: u.Quantity | None, **updates):
    """Clone cube metadata and preserve POLAROFF only when it is unambiguous."""
    if shared_polaroff is not None:
        updates["POLAROFF"] = shared_polaroff.to(u.degree)
    meta = clone_meta(cube, **updates)
    if shared_polaroff is None:
        meta.pop("POLAROFF", None)
    return meta


def _angle_for_scalar_meta(angle: u.Quantity) -> u.Quantity:
    """Return a scalar degree angle suitable for FITS header metadata."""
    angle = as_angle(angle, u.degree).to(u.degree)
    if angle.isscalar:
        return angle
    return np.nanmean(angle.to_value(u.degree)) * u.degree


def _warn_if_information_is_lost(pBp: np.ndarray, B: np.ndarray, context: str) -> None:
    """Warn when a lossy BP3 to BPB-style projection drops nonzero pBp."""
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


def _append_alpha(cubes: list[tuple[str, NDCube]], input_collection: NDCollection, mask: np.ndarray | None = None) -> list[tuple[str, NDCube]]:
    """Append alpha to an output cube list when the input collection had it."""
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


def _mzp_cubes_from_stack(
    data_stack: np.ndarray,
    input_collection: NDCollection,
    mask: np.ndarray | None = None,
    preferred_key: str | None = None,
) -> list[tuple[str, NDCube]]:
    """Convert a three-plane M/Z/P data stack into named NDCubes."""
    cube_template = get_template_cube(input_collection, preferred_key=preferred_key)
    shared_polaroff = _shared_polaroff(input_collection, get_data_keys(input_collection))
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


def _instrument_frame_polarizer_angles(input_collection: NDCollection) -> u.Quantity:
    """Return M/Z/P polarizer angles in the instrument frame."""
    data_shape = get_template_cube(input_collection, preferred_key="Z").data.shape
    return np.stack(
        [
            np.full(
                data_shape,
                wrap_linear_polarization(
                    as_angle(input_collection[key].meta["POLAR"], u.degree)
                    + as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree)
                ).to_value(u.degree),
            )
            * u.degree
            for key in ["M", "Z", "P"]
        ]
    )


def _solar_north_angles(input_collection: NDCollection) -> dict[str, u.Quantity]:
    """Return solar-north angles for each M/Z/P input cube."""
    data_shape = get_template_cube(input_collection, preferred_key="Z").data.shape
    lats = compute_lats(input_collection["Z"].wcs, data_shape)
    return {
        key: solnorth_from_wcs(input_collection[key].wcs, shape=data_shape, precomputed_lats=lats)
        for key in ["M", "Z", "P"]
    }


def _mzp_data_stack(collection: NDCollection) -> np.ndarray:
    """Return M/Z/P data in canonical order."""
    return np.stack([collection[key].data for key in ["M", "Z", "P"]], axis=0)


@transform(System.mzpsolar, System.bpb, use_alpha=True)
def mzpsolar_to_bpb(input_collection, **kwargs):
    """
    Notes
    ------
    Equation 7 and 9 in DeForest et al. 2022.

    """""
    alpha = _alpha_data(input_collection)
    polarizer_stack = np.stack(stack_data(input_collection, ["M", "Z", "P"]), axis=0)
    B, pB, pBp = bp3_from_polarizer_brightness(polarizer_stack, MZP_ANGLES, alpha)
    _warn_if_information_is_lost(pBp, B, "mzpsolar_to_bpb")

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="M")
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
    mzp_stack = bp3_to_polarizer_brightness(B, pB, np.zeros_like(pB), alpha, MZP_ANGLES)
    mask = combine_all_collection_masks(input_collection)
    cubes = _mzp_cubes_from_stack(mzp_stack, input_collection, mask=mask, preferred_key="B")
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bpb, System.btbr, use_alpha=True)
def bpb_to_btbr(input_collection, **kwargs):
    """Notes
    -----
    Equation 1 and 2 in DeForest et al. 2022.
    """
    B, pB = input_collection["B"].data, input_collection["pB"].data
    Br = (B - pB) / 2
    Bt = (B + pB) / 2

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="B")
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
    Bt, Br = input_collection["Bt"].data, input_collection["Br"].data
    pB = (Bt - Br)
    B = (Bt + Br)

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="Bt")
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

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="M")
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
        [1, -np.cos(2 * _angle_difference_radians(-60 * u.degree, alpha)), -np.sin(2 * _angle_difference_radians(-60 * u.degree, alpha))],
        [1, -np.cos(2 * _angle_difference_radians(0 * u.degree, alpha)), 0],
        [1, -np.cos(2 * _angle_difference_radians(60 * u.degree, alpha)), -np.sin(2 * _angle_difference_radians(60 * u.degree, alpha))],
    ])

    Bm = inv_mul_mx[0, 0] * Bi + inv_mul_mx[0, 1] * Bq + inv_mul_mx[0, 2] * Bu
    Bz = inv_mul_mx[1, 0] * Bi + inv_mul_mx[1, 1] * Bq + inv_mul_mx[1, 2] * Bu
    Bp = inv_mul_mx[2, 0] * Bi + inv_mul_mx[2, 1] * Bq + inv_mul_mx[2, 2] * Bu

    mask = combine_all_collection_masks(input_collection)
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
    polarizer_stack = np.stack(stack_data(input_collection, ["M", "Z", "P"]), axis=0)
    B, pB, pBp = bp3_from_polarizer_brightness(polarizer_stack, MZP_ANGLES, alpha)

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="M")
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
    mzp_stack = bp3_to_polarizer_brightness(B, pB, pBp, alpha, MZP_ANGLES)

    mask = combine_all_collection_masks(input_collection)
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
            Bt * np.sin(_angle_difference_radians(angle, alpha)) ** 2 + Br * np.cos(_angle_difference_radians(angle, alpha)) ** 2
            for angle in MZP_ANGLES
        ],
        axis=0,
    )
    mask = combine_all_collection_masks(input_collection)
    cubes = _mzp_cubes_from_stack(mzp_stack, input_collection, mask=mask, preferred_key="Bt")
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bp3, System.bthp, use_alpha=True)
def bp3_to_bthp(input_collection, **kwargs):
    """
    Notes
    ------
    Equations 9, 15, 16 in DeForest et al. 2022.

    Forward direction: (B, pB, pBp, alpha) -> (B, theta, p).

    theta is the polarization position angle measured CCW from solar north,
    wrapped to [-pi/2, pi/2) (linear-polarization axis convention).
    p is the fractional degree of linear polarization in [0, 1].
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

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="B")
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=template.meta)),
        ("theta", NDCube(theta, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Theta"))),
        ("p", NDCube(p, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Degree of Polarization"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
    return _collection_from_cubes(cubes)


@transform(System.bthp, System.bp3, use_alpha=True)
def bthp_to_bp3(input_collection, **kwargs):
    """
    Notes
    ------
    Inverse of Equations 15 and 16 in DeForest et al. 2022.

    Recovers (pB, pBp) from (theta, p, B, alpha) by reversing the
    bp3_to_bthp mapping::

        psi  = wrap_pm_pi(2 * (theta - pi/2 - alpha))
        pB   = p * B * cos(psi)
        pBp  = p * B * sin(psi)

    The wrapping to [-pi, pi) ensures the unique pre-image under
    wrap_linear_polarization used in the forward direction.
    """""
    B = input_collection["B"].data
    theta = as_angle(input_collection["theta"].data, u.radian).to_value(u.radian)
    p = input_collection["p"].data
    alpha = _alpha_data(input_collection).to_value(u.radian)

    psi = ((2 * (theta - np.pi / 2 - alpha)) + np.pi) % (2 * np.pi) - np.pi
    total_pol = p * B
    pB = total_pol * np.cos(psi)
    pBp = total_pol * np.sin(psi)

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="B")
    cubes = [
        ("B", NDCube(B, wcs=template.wcs, mask=mask, meta=template.meta)),
        ("pB", NDCube(pB, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="pB"))),
        ("pBp", NDCube(pBp, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="pB-prime"))),
    ]
    _append_alpha(cubes, input_collection, mask=mask)
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
    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="Bt")
    for angle in out_angles:
        value = Bt * np.sin(_angle_difference_radians(angle, alpha)) ** 2 + Br * np.cos(_angle_difference_radians(angle, alpha)) ** 2
        cubes.append((
            str(angle),
            NDCube(
                value,
                wcs=template.wcs,
                mask=mask,
                meta=clone_meta(template, POLAR=_angle_for_scalar_meta(angle)),
            ),
        ))
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
    input_keys = get_data_keys(input_collection)
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

    solved_stack = solve_to_three_polarizer_brightness(
        np.stack(stack_data(input_collection, input_keys), axis=0),
        observed_angles=phi,
        solved_angles=MZP_ANGLES,
        reference_angle=reference_angle,
    )
    mask = combine_all_collection_masks(input_collection)
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
    in_keys = get_data_keys(input_collection)
    source_angles = u.Quantity(
        [
            as_angle(input_collection[key].meta["POLAR"], u.degree)
            for key in in_keys
        ]
    )
    source_stack = np.stack(stack_data(input_collection, in_keys), axis=0)
    projected = solve_from_three_polarizer_brightness(
        source_stack,
        source_angles=source_angles,
        target_angles=out_angles,
        reference_angle=reference_angle,
    )

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection)
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
                    meta=_meta_with_shared_polaroff(
                        template,
                        shared_polaroff,
                        POLAR=_angle_for_scalar_meta(angle),
                    ),
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

    Sign convention (IAU / consistent with ``mzpsolar_to_stokes``)::

        I = B_0  + B_90
        Q = B_0  - B_90    (positive for N-S / along-north polarization)
        U = B_45 - B_135   (positive for NE-SW polarization)

    Q < 0 for tangential (east-west) Thomson-scattered coronal light,
    matching the sign produced by ``mzpsolar_to_stokes``.
    """""
    Bi = input_collection[str(0 * u.degree)].data + input_collection[str(90 * u.degree)].data
    Bq = input_collection[str(0 * u.degree)].data - input_collection[str(90 * u.degree)].data
    Bu = input_collection[str(45 * u.degree)].data - input_collection[str(135 * u.degree)].data

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key=str(0 * u.degree))
    cubes = [
        ("I", NDCube(Bi, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes I"))),
        ("Q", NDCube(Bq, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes Q"))),
        ("U", NDCube(Bu, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR="Stokes U"))),
    ]
    return _collection_from_cubes(cubes)


@transform(System.stokes, System.fourpol, use_alpha=False)
def stokes_to_fourpol(input_collection, **kwargs):
    """
    Notes
    ------
    Inverse of ``fourpol_to_stokes``.  Recovers four-polarizer brightness
    values from Stokes (I, Q, U) using the IAU sign convention::

        B_0   = (I + Q) / 2
        B_45  = (I + U) / 2
        B_90  = (I - Q) / 2
        B_135 = (I - U) / 2

    Roundtrip with ``fourpol_to_stokes`` is exact.
    """
    Bi = input_collection["I"].data
    Bq = input_collection["Q"].data
    Bu = input_collection["U"].data

    B0   = (Bi + Bq) / 2
    B45  = (Bi + Bu) / 2
    B90  = (Bi - Bq) / 2
    B135 = (Bi - Bu) / 2

    mask = combine_all_collection_masks(input_collection)
    template = get_template_cube(input_collection, preferred_key="I")
    pol_0   = str(0.0   * u.degree)
    pol_45  = str(45.0  * u.degree)
    pol_90  = str(90.0  * u.degree)
    pol_135 = str(135.0 * u.degree)
    cubes = [
        (pol_0,   NDCube(B0,   wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR=0.0   * u.degree))),
        (pol_45,  NDCube(B45,  wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR=45.0  * u.degree))),
        (pol_90,  NDCube(B90,  wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR=90.0  * u.degree))),
        (pol_135, NDCube(B135, wcs=template.wcs, mask=mask, meta=clone_meta(template, POLAR=135.0 * u.degree))),
    ]
    return _collection_from_cubes(cubes)


@transform(System.mzpsolar, System.mzpinstru, use_alpha=False)
def mzpsolar_to_mzpinstru(input_collection, reference_angle=0 * u.degree, **kwargs):
    """Convert solar-frame MZP to instrument-frame MZP.

    Notes
    -----
    Equation 45 in DeForest et al. 2022.

    **Angle convention**

    *Solar frame (mzpsolar)*
        - ``POLAR`` stores the solar position angle of each polarizer:
          M = −60°, Z = 0°, P = +60°, measured CCW from solar north.
        - ``POLARREF = "Solar"`` flags that ``POLAR`` is relative to solar north.
        - ``POLAROFF`` is an optional mechanical misalignment of the polarizer
          wheel *within* the instrument, not yet accounted for.

    *Instrument frame (mzpinstru)*
        - ``POLAR`` stores the same nominal solar-origin angle (−60 / 0 / 60).
        - ``POLARREF = "Instrument"`` flags that the *effective* angle in the
          image plane is ``POLAR + POLAROFF`` (i.e. after applying the wheel
          offset).
        - The brightness values themselves are what the physical polarizer
          (at instrument angle ``POLAR + POLAROFF``) would record, expressed
          in the coordinate system of the *instrument* rather than the Sun.

    **Projection step**
        For each target polarizer *j* (at instrument angle
        ``alpha_j = POLAR_j + POLAROFF_j``), build its solar-frame angle
        ``phi_j = alpha_j - solar_north`` and delegate the Equation 45
        projection to :func:`mzpsolar_to_npol`.
    """
    mask = combine_all_collection_masks(input_collection)
    solar_north = _solar_north_angles(input_collection)

    target_angles = []
    for key in ["M", "Z", "P"]:
        nominal_angle = as_angle(input_collection[key].meta["POLAR"], u.degree)
        polaroff = as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree)
        target_angles.append(wrap_linear_polarization(nominal_angle + polaroff - solar_north[key]))

    projected = mzpsolar_to_npol(
        input_collection,
        out_angles=u.Quantity(target_angles),
        reference_angle=reference_angle,
    )
    projected_stack = np.stack(stack_data(projected, get_data_keys(projected)), axis=0)

    cubes = []
    for key, value in zip(["M", "Z", "P"], projected_stack, strict=False):
        nominal_angle = as_angle(input_collection[key].meta["POLAR"], u.degree)
        polaroff = as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree)
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
    """Convert instrument-frame MZP to solar-frame MZP.

    Notes
    -----
    Equation 45 in DeForest et al. 2022.

    **Angle convention** — see :func:`mzpsolar_to_mzpinstru` for full
    definitions of ``POLAR``, ``POLAROFF``, and ``POLARREF``.

    **Projection step**
        The *effective instrument angle* of source polarizer *i* is::

            theta_i = POLAR_i + POLAROFF_i

        Convert those source angles into the solar frame with
        ``theta_i - solar_north`` and delegate the solve onto the canonical
        solar M/Z/P basis to :func:`npol_to_mzpsolar`.
    """
    mask = combine_all_collection_masks(input_collection)
    solar_north = _solar_north_angles(input_collection)

    input_keys = ["M", "Z", "P"]
    source_angles = u.Quantity(
        [
            wrap_linear_polarization(
                as_angle(input_collection[key].meta["POLAR"], u.degree)
                + as_angle(input_collection[key].meta.get("POLAROFF", 0), u.degree)
                - solar_north[key]
            )
            for key in input_keys
        ]
    )
    ordered_input = NDCollection(
        [(key, input_collection[key]) for key in input_keys],
        meta={},
        aligned_axes="all",
    )
    solved = npol_to_mzpsolar(
        ordered_input,
        in_angles=source_angles,
        reference_angle=reference_angle,
    )
    solved_stack = _mzp_data_stack(solved)

    cubes = []
    for key, mzp_angle, value in zip(["M", "Z", "P"], MZP_ANGLES, solved_stack, strict=False):
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
