"""Core transformation functions for solpolpy"""
import typing as t

import networkx as nx
from ndcube import NDCollection, NDCube

from solpolpy.alpha import radial_north
from solpolpy.constants import VALID_KINDS
from solpolpy.errors import UnsupportedTransformationError
from solpolpy.graph import transform_graph
from solpolpy.instruments import load_data
from solpolpy.polarizers import npol_to_mzp


def resolve(input_data: t.Union[t.List[str], NDCollection], out_system: str) -> NDCollection:
    """
    Apply - apply a polarization transformation to a set of input
    dataframes.

    Parameters
    ----------
    input_data : NDCollection or List[str]
        Either: 1) a collection where each member NDCube has an expected name or 2) a list of paths to FITS files.
        We recommend option 2.

    out_system : string
        The polarization state you want to convert your input dataframes to.
        Must be one of the following strings:

        - "MZP": Triplet of images taken at -60°, 0°, and +60° polarizing angles.
        - "BtBr": Pair of images with polarization along the tangential and radial direction with respect to the Sun respectively.
        - "Stokes": Total brightness ("I"), polarized brightness along vertical and horizontal axes (Q) and polarized brightness along ±45° (U) .
        - "BpB": Total brightness and ‘excess polarized’ brightness images pair respectively.
        - "Bp3": Analogous to Stokes I, Q and U, but rotates around the Sun instead of a fixed frame of reference of the instrument.
        - "Bthp": Total brightness, angle and degree of polarization.
        - "fourpol": For observations taken at sequence of four polarizer angles, i.e. 0°, 45°, 90° and 135°.
        - "npol": Set of images taken at than three polarizing angles other than MZP

    Raises
    ------
    AssertionError
      This gets raised if the data cannot be converted or polarization
      transformation cannot be calculated due to a discontinuity or infinity.

    Returns
    -------
    NDCollection
        The transformed data are returned as a NDcollection.  Most
        Transforms maintain the dimensionality of the source vectors.  Some
        embed (increase dimensionality of the vectors) or project (decrease
        dimensionality of the vectors); additional input dimensions, if
        present, are still appended to the output vectors in all any case.
    """
    if isinstance(input_data, list):
        input_data = load_data(input_data)

    input_kind = determine_input_kind(input_data)

    # if it's npol we immediately standardize to MZP
    if input_kind == "npol":
        input_data = npol_to_mzp(input_data)
        input_kind = "MZP"

    input_key = list(input_data)
    transform_path = get_transform_path(input_kind, out_system)
    equation = get_transform_equation(transform_path)
    requires_alpha = check_alpha_requirement(transform_path)

    if requires_alpha and "alpha" not in input_key:
        input_data = add_alpha(input_data)
    result = equation(input_data)

    return result


def determine_input_kind(input_data: NDCollection) -> str:
    """Determine what kind of data was input in the NDCollection

    Parameters
    ----------
    input_data : NDCollection
        data to evaluate kind of

    Returns
    -------
    str
        a valid input kind, see documentation of `resolve` for the full list under `out_system`
    """
    input_keys = list(input_data)
    for valid_kind, param_list in VALID_KINDS.items():
        for param_option in param_list:
            if set(input_keys) == set(param_option):
                return valid_kind
    raise ValueError("Unidentified Polarization System.")


def get_transform_path(input_kind: str, output_kind: str) -> t.List[str]:
    """Given an input and output system type, determine the require path of transforms from the transform graph

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
        raise UnsupportedTransformationError(f"Not possible to convert {input_kind} to {output_kind}")
    return path


def get_transform_equation(path: t.List[str]) -> t.Callable:
    """Given a transform path, compose the equation, i.e. the composed function of transforms that executes that path

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
        current_function = _compose2(transform_graph.get_edge_data(step_start, step_end)['func'],
                                     current_function)
    return current_function


def check_alpha_requirement(path: t.List[str]) -> bool:
    """Determine if an alpha array is required for this transformation path

    Parameters
    ----------
    path : List[str]
        a path through the transform graph

    Returns
    -------
    bool
        whether alpha array is required for transformation
    """
    requires_alpha = False
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i + 1]
        requires_alpha = transform_graph.get_edge_data(step_start, step_end)['requires_alpha'] or requires_alpha
    return requires_alpha


def _determine_image_shape(input_collection: NDCollection) -> t.Tuple[int, int]:
    """Evaluates the shape of the image in the input NDCollection

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
    shape = input_collection[keys[0]].data.shape
    return shape


def add_alpha(input_data: NDCollection) -> NDCollection:
    """Adds an alpha array to an image NDCollection

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
    keys = list(input_data)
    wcs = input_data[keys[0]].wcs
    metad = input_data[keys[0]].meta

    if len(img_shape) == 2:  # it's an image and not just an array
        alpha = radial_north(img_shape)
    else:
        raise ValueError(f"Data must be an image with 2 dimensions, found {len(img_shape)}.")
    input_data.update(NDCollection([("alpha", NDCube(alpha, wcs=wcs, meta=metad))], meta={}, aligned_axes='all'))

    return input_data


def _compose2(f: t.Callable, g: t.Callable) -> t.Callable:
    """Compose 2 functions together, i.e. f(g(x))

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
    return lambda *a, **kw: f(g(*a, **kw))


def identity(x: t.Any) -> t.Any:
    """Identity function that returns the input

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
