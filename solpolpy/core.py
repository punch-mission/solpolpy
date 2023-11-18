from typing import Callable, List, Tuple, Union
import numbers

from astropy.io import fits
import astropy.units as u
from astropy.units.quantity import Quantity
import numpy as np
import networkx as nx
import warnings
import solpolpy as sp

from ndcube import NDCube, NDCollection
from astropy.wcs import WCS

from solpolpy.constants import VALID_KINDS
from solpolpy.graph import transform_graph
from solpolpy.alpha import radial_north, radial_west
from solpolpy.instruments import load_data


def resolve(input_data: Union[List[str], NDCollection], out_polarize_state: str) -> NDCollection:
    """
    Apply - apply a polarization transformation to a set of input
    dataframes.

    Parameters
    ----------
    input_data : NDCollection
        NDCollection formatted as follows:

            - Stokes NDCollection
                "Bi":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "Bq":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "Bu":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "Bv":np.array - Should be included as a triplet of I,Q,U, and optionally V

            - Brightness & Polarized brightness NDCollection
                "B":np.array - Should be included as a double of B, pB
                "pB":np.array - Should be included as a doublet of B, pB

            - Radial & tangential brightness NDCollection
                "Br":np.array - Should be included as a doublet of Bt, Br
                "Bt":np.array - Should be included as a doublet of Bt, Br

            - MZP triplet NDCollection
                "Bm":np.array - Should be included as a triplet of M,Z,P
                "Bz":np.array - Should be included as a triplet of M,Z,P
                "Bp":np.array - Should be included as a triplet of M,Z,P

            - Angular (npol) NDCollection [where keys are ints or floats representing
            polarizer angle]
                X1:np.array[xn] - Should be included as at least a triplet of angles
                X2:np.array[xn] - Should be included as at least a triplet of angles
                X3:np.array[xn] - Should be included as at least a triplet of angles
                ...                               ...
                Xn:np.array[xn] - Should be included as at least a triplet of angles

    out_polarize_state : string
      This is the polarization state you want to convert your input
      dataframes to.  These include:

        - stokes: convert the input dataframes to output Stokes dataframes. If
            3 polarization angles are input, or can be derived, the I, Q, and U
            Stokes parameters are output. If 4 or more polarization angles are
            provided, and the angles are conducsive, the full I,Q,U, and V
            Stokes parameters are provided.

        - B: converts input dataframes into B ("unpolarized
            brightness") parameters.

        - pB: converts input dataframes into pB ("polarized brightness")
            parameters.

        - Bt: Produces "Bt", the radiance observed through a linear polarizer
            oriented tangentially to a solar-concentric circle passing through
            the image point of interest.

        - Br: Produces "Br" the radiance observed through a linear polarizer
            oriented radially to the centre Sun through the same point.

        - BrBt: Produces both the  radiance observed through a linear polarizer
            oriented tangentially to a solar-concentric circle passing through
            the image point of interest "Bt", and the radiance observed through
            a linear polarizer oriented radially to the centre Sun through the
            same point "Br".

        - MZP: converts input dataframes into a system of virtual
            polarizer triplets each separated by 60 degrees (Minus [-60],
            Zero [0], Plus [60]).

        - 4pol: converts input dataframes into a system of virtual polarizer
            triplets each separated by 45 degrees at -45, 0, 45, 90.

        - Xpol: converts input dataframes into a system of virtual polarizer
            triplets each separated by X degrees, specified by the optional
            input separation, which should be input in degrees.

        - BpB: converts input dataframes into two dataframes, B ("unpolarized
            brightness") parameter and its counterpart pB ("polarized brightness")

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
    deg2rad = (np.pi * u.radian) / (180 * u.degree)

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
