from typing import Callable, List, Tuple, Union

import astropy.units as u
import numpy as np
import networkx as nx
import solpolpy as sp
from ndcube import NDCube, NDCollection

from solpolpy.constants import VALID_KINDS
from solpolpy.graph import transform_graph
from solpolpy.alpha import radial_north
from solpolpy.instruments import load_data


def resolve(input_data: Union[List[str], NDCollection], out_polarize_state: str) -> NDCollection:
    """
    Apply - apply a polarization transformation to a set of input
    dataframes.

    Parameters
    ----------
    input_data : NDCollection or List[str]
        Either: 1) a collection where each member NDCube has an expected name or 2) a list of paths to FITS files.
        We recommend option 2.

    out_polarize_state : string
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
    if input_kind == "npol":
        input_data = sp.polarizers.npol_to_mzp(input_data)
        input_kind = "MZP"

    input_key = list(input_data)
    transform_path = get_transform_path(input_kind, out_polarize_state)
    equation = get_transform_equation(transform_path)
    requires_alpha = check_alpha_requirement(transform_path)

    if requires_alpha and "alpha" not in input_key:
        input_data = add_alpha(input_data)
    result = equation(input_data)

    return result


def determine_input_kind(input_data: NDCollection) -> str:
    input_keys = list(input_data)
    for valid_kind, param_list in VALID_KINDS.items():
        for param_option in param_list:
            if set(input_keys) == set(param_option):
                return valid_kind
    raise ValueError("Unidentified Polarization System.")


def get_transform_path(input_kind: str, output_kind: str) -> List[str]:
    try:
        path = nx.shortest_path(transform_graph, input_kind, output_kind)
    except nx.exception.NetworkXNoPath:
        raise RuntimeError(f"Not possible to convert {input_kind} to {output_kind}")  # TODO: make this a custom error
    return path


def get_transform_equation(path: List[str]) -> Callable:
    current_function = identity
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i + 1]
        current_function = _compose2(transform_graph.get_edge_data(step_start, step_end)['func'],
                                     current_function)
    return current_function


def check_alpha_requirement(path: List[str]):
    requires_alpha = False
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i + 1]
        requires_alpha = transform_graph.get_edge_data(step_start, step_end)['requires_alpha'] or requires_alpha
    return requires_alpha


def _determine_image_shape(input_cube) -> Tuple[int, int]:
    keys = list(input_cube.keys())
    shape = input_cube[keys[0]].data.shape
    return shape


def add_alpha(input_data: NDCollection) -> NDCollection:
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


def _compose2(f, g):
    return lambda *a, **kw: f(g(*a, **kw))


def identity(x):
    return x
