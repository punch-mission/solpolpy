from typing import Callable, Dict, List, Tuple
import functools
from collections import Counter

from astropy.io import fits
import numpy as np
import networkx as nx

from solpolpy.constants import VALID_KINDS
from solpolpy.graph import transform_graph
from solpolpy.alpha import ALPHA_FUNCTIONS


def resolve(input_data, out_polarize_state, separation=None, alpha=None, Error=False):
    """
    Apply - apply a depolarization transform to a set of input
    dataframes.

    Parameters
    ----------
    input_data : Dict[str, np.ndarray]
        Dictionary formatted as follows:

            - Stokes dictionary
                "I":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "Q":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "U":np.array - Should be included as a triplet of I,Q,U, and optionally V
                "V":np.array - Should be included as a triplet of I,Q,U, and optionally V

            - Brightness & Polarized brightness dictionary
                "B":np.array - Should be included as a double of B, pB
                "pB":np.array - Should be included as a doublet of B, pB

            - Radial & tangential brightness dictionary
                "Br":np.array - Should be included as a doublet of Bt, Br
                "Bt":np.array - Should be included as a doublet of Bt, Br

            - MZP triplet dictionary
                "M":np.array - Should be included as a triplet of M,Z,P
                "Z":np.array - Should be included as a triplet of M,Z,P
                "P":np.array - Should be included as a triplet of M,Z,P

            - Angular dictionary [where keys are ints or floats representing
            polarizer angle]
                X1:np.array[xn] - Should be included as at least a triplet of angles
                X2:np.array[xn] - Should be included as at least a triplet of angles
                X3:np.array[xn] - Should be included as at least a triplet of angles
                ...                               ...
                Xn:np.array[xn] - Should be included as at least a triplet of angles

    out_polarize_state : string
      This is the polarization state you want to convert your input
      dataframes to.  These include:

        - Stokes: convert the input dataframes to output Stokes dataframes. If
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

        - 3pol or MZP: converts input dataframes into a system of virtual
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
    numpy.ndarray
        The transformed vector data are returned as a numpy.ndarray.  Most
        Transforms maintain the dimensionality of the source vectors.  Some
        embed (increase dimensionality of the vectors) or project (decrease
        dimensionality of the vectors); additional input dimensions, if
        present, are still appended to the output vectors in all any case.
    """
    if isinstance(input_data, list):
        input_data = convert_image_list_to_dict(input_data)
    input_kind = determine_input_kind(input_data)
    input_data = _add_alpha(input_data, alpha)
    equation = get_transform_equation(input_kind, out_polarize_state)
    result = equation(input_data)
    # uses_alpha = get_alpha_usage(input_kind, out_polarize_state)
    # if uses_alpha:
    #     result = equation(input_data, alpha)
    # else:
    #     result = equation(input_data)
    return result


def determine_input_kind(input_data: Dict[str, np.ndarray]) -> str:
    input_keys = list(input_data.keys())
    for valid_kind, param_list in VALID_KINDS.items():
        for param_option in param_list:
            if set(input_keys) == set(param_option):
                return valid_kind
    input_key_kinds = Counter([type(v) for v in input_keys])
    if input_key_kinds[np.float64] >= 3:
        return "MZP"
    raise ValueError("Invalid Keys")  # TODO: improve this message


def get_transform_equation(input_kind: str, output_kind: str) -> Callable:
    try:
        path = nx.shortest_path(transform_graph, input_kind, output_kind)
    except nx.exception.NetworkXNoPath:
        raise RuntimeError(f"Not possible to convert {input_kind} to {output_kind}") # TODO: make this a custom error

    current_function = identity
    for i, step_start in enumerate(path):
        if i+1 < len(path):
            step_end = path[i+1]
            current_function = _compose2(transform_graph.get_edge_data(step_start, step_end)['func'],
                                         current_function)
    return current_function


def _determine_image_shape(input_dict: Dict[str, np.ndarray]) -> Tuple[int, int]:
    keys = list(input_dict.keys())
    shape = input_dict[keys[0]].shape
    return shape


def _add_alpha(input_data: Dict[str, np.ndarray], kind) -> Dict[str, np.ndarray]:
    # test if alpha exists. if not check if alpha keyword added. if not create default alpha with warning.

    img_shape = _determine_image_shape(input_data)
    if len(img_shape) == 2:  # it's an image and not just an array
        kind = "radial" if kind is None else kind
        if "alpha" not in input_data:
            if kind not in ALPHA_FUNCTIONS:
                raise ValueError(f"Requested a {kind} alpha type. This is not valid. Must be in {ALPHA_FUNCTIONS.keys()}")
            input_data['alpha'] = ALPHA_FUNCTIONS[kind](img_shape)
    return input_data


def _convert_STEREO_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    # data_out={}
    #     list_len=len(data_in)
    #     assert list_len >= 2, 'requires at least 2 FITS files'
    #
    #     for xlist_item in data_in:
    #         with fits.open(xlist_item) as hdul:
    #             assert hdul[0].header['INSTRUME'] == 'SECCHI', 'requires FITS to be SECCHI COR data files'
    #             data_out[hdul[0].header['POLAR']]=hdul[0].data
    #             image_hdr = hdul[0].header
    pass


def _convert_LASCO_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    # data_out={}
    #     list_len=len(data_in)
    #     assert list_len >= 2, 'requires at least 2 FITS files'
    #
    #     for xlist_item in data_in:
    #         with fits.open(xlist_item) as hdul:
    #             assert hdul[0].header['INSTRUME'] == 'SECCHI', 'requires FITS to be SECCHI COR data files'
    #             data_out[hdul[0].header['POLAR']]=hdul[0].data
    #             image_hdr = hdul[0].header
    pass


def convert_image_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    # data_out={}
    #     list_len=len(data_in)
    #     assert list_len >= 2, 'requires at least 2 FITS files'
    #
    #     for xlist_item in data_in:
    #         with fits.open(xlist_item) as hdul:
    #             assert hdul[0].header['INSTRUME'] == 'SECCHI', 'requires FITS to be SECCHI COR data files'
    #             data_out[hdul[0].header['POLAR']]=hdul[0].data
    #             image_hdr = hdul[0].header
    # if input_data is a STEREO:
    #     out = _convert_STEREO_list_to_dict(input_data)
    # elif input_data is LASCO:
    #     out = _convert_LASCO_list_to_dict(input_data)
    # else:
    #     raise Exception("Don't recognize this FITS type. Use dictionary input.")
    pass


def get_alpha_usage(input_kind: str, output_kind: str) -> bool:
    pass


def _compose2(f, g):
    return lambda *a, **kw: f(g(*a, **kw))


def identity(x):
    return x
