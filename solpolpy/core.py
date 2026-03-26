"""Core transformation functions for solpolpy."""
import typing as t
import warnings

import astropy.units as u
import networkx as nx
from ndcube import NDCollection, NDCube

from solpolpy.alpha import radial_from_wcs, radial_north
from solpolpy.constants import STEREOA_REFERENCE_ANGLE, STEREOB_REFERENCE_ANGLE
from solpolpy.errors import UnsupportedTransformationError
from solpolpy.instruments import load_data
from solpolpy.physics import wrap_pm_pi
from solpolpy.transforms import SYSTEM_REQUIRED_KEYS, System, transform_graph
from solpolpy.util import solnorth_from_wcs


@u.quantity_input
def resolve(input_data: list[str] | NDCollection,
            out_system: str,
            in_angles: u.degree = None,
            out_angles: u.degree = None,
            reference_angle: u.degree = None) -> NDCollection:
    """Apply a polarization transformation to a set of input dataframes.

    Parameters
    ----------
    input_data : NDCollection or List[str]
        Either: 1) a collection where each member NDCube has an expected name or 2) a list of paths to FITS files.
        We recommend option 2.

    out_system : string
        The polarization state you want to convert your input dataframes to.
        Must be one of the following strings:

        - "mzpsolar": Triplet of images taken at -60°, 0°, and +60° polarizing angles with a reference angle set to solar frame.
        - "mzpinstru": Triplet of images taken at -60°, 0°, and +60° polarizing angles with a reference angle set to instrument frame.
        - "btbr": A Pair of images with polarization along the tangential and radial direction with respect to the Sun respectively.
        - "stokes": Total brightness ("I"), polarized brightness along vertical and horizontal axes (Q) and polarized brightness along ±45° (U) .
        - "bpb": Total brightness and ‘excess polarized’ brightness images pair respectively.
        - "bp3": Analogous to Stokes I, Q and U, but rotates around the Sun instead of a fixed frame of reference of the instrument.
        - "bthp": Total brightness, angle and degree of polarization.
        - "fourpol": For observations taken at sequence of four polarizer angles, i.e. 0°, 45°, 90° and 135°.
        - "npol": Set of images taken at than arbitrary polarizing angles other than MZP

    out_angles : u.degree
        Angles to use when converting to npol or some arbitrary system

    in_angles : u.degree
        Angles to use when converting to from npol or some arbitrary system to mzpsolar

    reference_angle : u.degree
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

    if getattr(equation, "requires_out_angles", False) and out_angles is None:
        raise ValueError("Out angles must be specified for this transform.")

    if getattr(equation, "requires_in_angles", False) and in_angles is None:
        raise ValueError("In angles must be specified for this transform.")

    if requires_alpha(equation) and "alpha" not in input_keys:
        input_data = add_alpha(input_data)

    reference_angle = determine_reference_angle(input_data) if reference_angle is None else reference_angle

    return equation(input_data,
                    reference_angle=reference_angle,
                    in_angles=in_angles,
                    out_angles=out_angles)


def determine_reference_angle(input_collection: NDCollection) -> u.degree:
    """Get the instrument specific offset angle."""
    first_key = next(iter(input_collection.keys()))
    match input_collection[first_key].meta.get("OBSRVTRY", "BLANK"):
        case "STEREO_A":
            reference_angle = STEREOA_REFERENCE_ANGLE
        case "STEREO_B":
            reference_angle = STEREOB_REFERENCE_ANGLE
        case _:
            reference_angle = 0 * u.degree

    return reference_angle


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
        if valid_kind in [System.mzpinstru, System.mzpsolar] and param_set == input_keys:
            polarref_value = input_data['Z'].meta.get("POLARREF", "solar").lower()
            return System.mzpinstru if polarref_value == "instrument" else System.mzpsolar
        if valid_kind != System.npol and param_set == input_keys:
            return valid_kind
    try:
        input_keys_quantities = [u.Quantity(key) for key in input_keys]
    except (TypeError, ValueError):
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
        try:
            alpha = radial_from_wcs(wcs, img_shape)
        except Exception as err:  # pragma: no cover - best-effort fallback
            alpha = radial_north(img_shape)
            try:
                alpha = wrap_pm_pi(alpha + solnorth_from_wcs(wcs, img_shape).to(u.radian))
            except Exception:
                pass
            warnings.warn(
                f"Falling back to image-centered alpha approximation because WCS-based solar-center calculation failed: {err}",
                stacklevel=2,
            )
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
    out.uses_in_angles = getattr(f, "uses_in_angles", False) or getattr(g, "uses_in_angles", False)
    out.requires_out_angles = getattr(f, "requires_out_angles", False) or getattr(g, "requires_out_angles", False)
    out.requires_in_angles = getattr(f, "requires_in_angles", False) or getattr(g, "requires_in_angles", False)

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
