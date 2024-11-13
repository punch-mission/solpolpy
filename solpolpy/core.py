"""Core transformation functions for solpolpy."""
import typing as t

import astropy.units as u
import networkx as nx
import numpy as np
from ndcube import NDCollection, NDCube

from solpolpy.alpha import radial_north
from solpolpy.constants import STEREOA_OFFSET_ANGLE, STEREOB_OFFSET_ANGLE
from solpolpy.errors import UnsupportedTransformationError
from solpolpy.instruments import load_data
from solpolpy.transforms import SYSTEM_REQUIRED_KEYS, System, transform_graph


@u.quantity_input
def resolve(input_data: list[str] | NDCollection,
            out_system: str,
            imax_effect: bool = False,
            out_angles: u.degree = None,
            offset_angle: u.degree = None) -> NDCollection:
    """Apply a polarization transformation to a set of input dataframes.

    Parameters
    ----------
    input_data : NDCollection or List[str]
        Either: 1) a collection where each member NDCube has an expected name or 2) a list of paths to FITS files.
        We recommend option 2.

    out_system : string
        The polarization state you want to convert your input dataframes to.
        Must be one of the following strings:

        - "mzp": Triplet of images taken at -60°, 0°, and +60° polarizing angles.
        - "btbr": Pair of images with polarization along the tangential and radial direction with respect to the Sun respectively.
        - "stokes": Total brightness ("I"), polarized brightness along vertical and horizontal axes (Q) and polarized brightness along ±45° (U) .
        - "bpb": Total brightness and ‘excess polarized’ brightness images pair respectively.
        - "bp3": Analogous to Stokes I, Q and U, but rotates around the Sun instead of a fixed frame of reference of the instrument.
        - "bthp": Total brightness, angle and degree of polarization.
        - "fourpol": For observations taken at sequence of four polarizer angles, i.e. 0°, 45°, 90° and 135°.
        - "npol": Set of images taken at than arbitrary polarizing angles other than MZP

    imax_effect : Boolean
        The 'IMAX effect' describes the change in apparent measured polarization angle as an result of foreshortening effects.
        This effect becomes more pronounced for wide field polarized imagers - see Patel et al (2024, in preparation)
        If True, applies the IMAX effect for wide field imagers as part of the resolution process.

    out_angles : u.degree
        Angles to use when converting to npol or some arbitrary system

    offset_angle : u.degree
        Reference angle used for the polarizer offset. If None, it will try to determine it from the metadata.

    Returns
    -------
    NDCollection
        The transformed data are returned as a NDCollection.

    """
    out_system = out_system.lower()

    if isinstance(input_data, list):
        input_data = load_data(input_data)

    input_kind = determine_input_kind(input_data)

    input_keys = list(input_data.keys())
    transform_path = get_transform_path(input_kind, out_system)
    equation = get_transform_equation(transform_path)

    if getattr(equation, "uses_out_angles", False) and out_angles is None:
        raise ValueError("Out angles must be specified for this transform.")

    if imax_effect:
        if input_kind == "mzp":
            input_data = resolve_imax_effect(input_data)
        else:
            msg = "IMAX effect applies only for transformations starting from MZP."
            raise UnsupportedTransformationError(msg)

    if requires_alpha(equation) and "alpha" not in input_keys:
        input_data = add_alpha(input_data)

    offset_angle = determine_offset_angle(input_data) if offset_angle is None else offset_angle

    return equation(input_data,
                    offset_angle=offset_angle,
                    out_angles=out_angles)


def determine_offset_angle(input_collection: NDCollection) -> u.degree:
    """Get the instrument specific offset angle."""
    first_key = next(iter(input_collection.keys()))
    match input_collection[first_key].meta.get("OBSRVTRY", "BLANK"):
        case "STEREO_A":
            offset_angle = STEREOA_OFFSET_ANGLE
        case "STEREO_B":
            offset_angle = STEREOB_OFFSET_ANGLE
        case _:
            offset_angle = 0 * u.degree

    return offset_angle


def determine_input_kind(input_data: NDCollection) -> System:
    """Determine what kind of data was input in the NDCollection.

    Parameters
    ----------
    input_data : NDCollection
        data to evaluate kind of

    Returns
    -------
    str
        a valid input kind, see documentation of `resolve` for the full list under `out_system`

    """
    input_keys = set(input_data)
    input_keys.discard("alpha")
    if len(input_keys) == 0:
        msg = "Found no cubes in the `input_data` collection."
        raise ValueError(msg)

    for valid_kind, param_set in SYSTEM_REQUIRED_KEYS.items():
        if valid_kind != System.npol and param_set == input_keys:
            return valid_kind
    try:
        input_keys_quantities = [u.Quantity(key) for key in input_keys]
    except TypeError:
        pass
    else:
        if all(u.get_physical_type(q) == "angle" for q in input_keys_quantities):
            return System.npol

    msg = "Could not determine input transformation."
    raise UnsupportedTransformationError(msg)


def get_transform_path(input_kind: str, output_kind: str) -> list[str]:
    """Given an input and output system type, determine the require path of transforms from the transform graph.

    Parameters
    ----------
    input_kind : str
        starting point for transformations

    output_kind : str
        ending point for transformations

    Returns
    -------
    List[str]
        a list of transformation identifiers used to convert from `input_kind` to `output_kind`

    """
    try:
        path = nx.shortest_path(transform_graph, input_kind, output_kind)
    except nx.exception.NetworkXNoPath:
        msg = f"Not possible to convert {input_kind} to {output_kind}"
        raise UnsupportedTransformationError(msg)
    return path


def get_transform_equation(path: list[str]) -> t.Callable:
    """Given a transform path, compose the equation, i.e. the composed function of transforms that executes that path.

    Parameters
    ----------
    path : List[str]
        a list of transform identifiers from the path

    Returns
    -------
    Callable
        a function that executes the transformation

    """
    current_function = identity
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i + 1]
        current_function = _compose2(transform_graph.get_edge_data(step_start, step_end)["func"],
                                     current_function)
    return current_function


def requires_alpha(func: t.Callable) -> bool:
    """Determine if an alpha array is required for this transformation path.

    Parameters
    ----------
    path : List[str]
        a path through the transform graph

    Returns
    -------
    bool
        whether alpha array is required for transformation

    """
    return getattr(func, "uses_alpha", False)


def generate_imax_matrix(array_shape) -> np.ndarray:
    """Define an A matrix with which to convert MZP^ (camera coords) = A x MZP (solar coords).

    Parameters
    ----------
    array_shape
        Defined input WCS array shape for matrix generation

    Returns
    -------
    ndarray
        Output A matrix used in converting between camera coordinates and solar coordinates

    """
    # Ideal MZP wrt Solar North
    # TODO fix variable name
    thmzp = [-60, 0, 60] * u.degree

    # TODO don't hard code extents but get from the file
    long_arr, lat_arr = np.meshgrid(np.linspace(-20, 20, array_shape[0]), np.linspace(-20, 20, array_shape[1]))

    # Foreshortening (IMAX) effect on polarizer angle
    phi_m = np.arctan2(np.tan(thmzp[0]) * np.cos(long_arr * u.degree), np.cos(lat_arr * u.degree)).to(u.degree)
    phi_z = np.arctan2(np.tan(thmzp[1]) * np.cos(long_arr * u.degree), np.cos(lat_arr * u.degree)).to(u.degree)
    phi_p = np.arctan2(np.tan(thmzp[2]) * np.cos(long_arr * u.degree), np.cos(lat_arr * u.degree)).to(u.degree)

    phi = np.stack([phi_m, phi_z, phi_p])

    # Define the A matrix
    mat_a = np.empty((array_shape[0], array_shape[1], 3, 3))

    for i in range(3):
        for j in range(3):
            mat_a[:, :, i, j] = (4 * np.cos(phi[i] - thmzp[j]) ** 2 - 1) / 3

    return mat_a


def resolve_imax_effect(input_data: NDCollection) -> NDCollection:
    """Resolves the IMAX effect for provided input data, correcting measured polarization angles for wide FOV imagers.

    Parameters
    ----------
    input_data : NDCollection
        Input data on which to correct foreshortened polarization angles

    Returns
    -------
    NDCollection
        Output data with corrected polarization angles

    """
    data_shape = _determine_image_shape(input_data)
    data_mzp_camera = np.zeros([data_shape[0], data_shape[1], 3, 1])

    for i, key in enumerate(["M", "Z", "P"]):
        data_mzp_camera[:, :, i, 0] = input_data[key].data

    imax_matrix = generate_imax_matrix(data_shape)
    try:
        imax_matrix_inv = np.linalg.inv(imax_matrix)
    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            msg = "Singular IMAX effect matrix is degenerate"
            raise ValueError(msg)
        else:
            raise

    data_mzp_solar = np.matmul(imax_matrix_inv, data_mzp_camera)

    for i, key in enumerate(input_data.keys()):
        input_data[key].data[:, :] = data_mzp_solar[:, :, i, 0]

    return input_data


def _determine_image_shape(input_collection: NDCollection) -> tuple[int, int]:
    """Evaluates the shape of the image in the input NDCollection.

    Parameters
    ----------
    input_collection : NDCollection
        collection to determine image shape for

    Returns
    -------
    Tuple[int, int]
        shape of image data

    """
    keys = list(input_collection.keys())
    return input_collection[keys[0]].data.shape


def add_alpha(input_data: NDCollection) -> NDCollection:
    """Adds an alpha array to an image NDCollection.

    Parameters
    ----------
    input_data : NDCollection
        dataset to append alpha to

    Returns
    -------
    NDCollection
        dataset with alpha array appended

    """
    # test if alpha exists. if not check if alpha keyword added. if not create default alpha with warning.
    img_shape = _determine_image_shape(input_data)
    keys = list(input_data.keys())
    wcs = input_data[keys[0]].wcs
    meta = input_data[keys[0]].meta

    if len(img_shape) == 2:  # it's an image and not just an array
        alpha = radial_north(img_shape)
    else:
        msg = f"Data must be an image with 2 dimensions, found {len(img_shape)}."
        raise ValueError(msg)
    input_data.update(NDCollection([("alpha", NDCube(alpha, wcs=wcs, meta=meta))], meta={},
                                   aligned_axes="all"))

    return input_data


def _compose2(f: t.Callable, g: t.Callable) -> t.Callable:
    """Compose 2 functions together, i.e. f(g(x)).

    Parameters
    ----------
    f : Callable
        outer function
    g : Callable
        inner function

    Returns
    -------
    Callable
        composed function

    """
    def out(*a, **kw):
        return f(g(*a, **kw), **kw)
    out.uses_alpha = getattr(f, "uses_alpha", False) or getattr(g, "uses_alpha", False)
    out.uses_out_angles = getattr(f, "uses_out_angles", False) or getattr(g, "uses_out_angles", False)

    return out


def identity(x: t.Any, **kwargs) -> t.Any:
    """Identity function that returns the input.

    Parameters
    ----------
    x : Any
        value

    Returns
    -------
    Any
        input value returned back

    """
    return x
