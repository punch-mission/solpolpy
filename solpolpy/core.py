from typing import Callable, Dict, List, Tuple
import functools
from collections import Counter
import numbers

from astropy.io import fits
import astropy.units as u
from astropy.units.quantity import Quantity
import numpy as np
import networkx as nx

from solpolpy.constants import VALID_KINDS
from solpolpy.graph import transform_graph
from solpolpy.alpha import ALPHA_FUNCTIONS


def resolve(input_data, out_polarize_state, alpha=None):
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
    input_data, input_has_radians = sanitize_data_dict(input_data, u.radian)
    input_kind = determine_input_kind(input_data)

    transform_path = get_transform_path(input_kind, out_polarize_state)
    equation = get_transform_equation(transform_path)
    requires_alpha = check_alpha_requirement(transform_path)
    if requires_alpha:
        input_data = add_alpha(input_data, alpha)

    result = equation(input_data)
    if not input_has_radians:
        result, _ = sanitize_data_dict(result, u.degree)
    return result


def sanitize_data_dict(input_data, output_type):
    if output_type not in [u.radian, u.degree]:
        raise RuntimeError(f"Keys must be converted to degrees or radians. Found ouput_type={output_type}")

    # Set up necessary conversions and
    deg2rad = (np.pi * u.radian) / (180 * u.degree)
    rad2deg = (180 * u.degree) / (np.pi * u.radian)

    # Sanitize all the keys
    found_radians = False
    output_dict = dict()
    for key, value in input_data.items():
        if isinstance(key, numbers.Real):
            if output_type == u.degree:
                output_dict[key * u.degree] = value
            else:
                output_dict[np.deg2rad(key) * u.radian] = value
        elif isinstance(key, Quantity):
            if key.unit not in [u.radian, u.degree]:
                raise RuntimeError(f"Unsupported key type of {u.key} found. Must be radians or degrees.")
            if key.unit == u.radian:
                found_radians = True
                if output_type == u.degree:
                    output_dict[key*rad2deg] = value
                else:
                    output_dict[key] = value
            elif key.unit == u.degree:
                if output_type == u.radian:
                    output_dict[key*deg2rad] = value
                else:
                    output_dict[key] = value
        else:
            output_dict[key] = value

    # Sanitize alpha map
    if "alpha" in input_data:
        if isinstance(input_data['alpha'], Quantity):
            if input_data['alpha'].unit == u.radian:
                if output_type == u.degree:
                    output_dict["alpha"] = input_data['alpha'] * rad2deg
                else:
                    output_dict["alpha"] = input_data["alpha"]
            elif input_data['alpha'].unit == u.degree:
                if output_type == u.degree:
                    output_dict["alpha"] = input_data["alpha"]
                else:
                    output_dict["alpha"] = input_data["alpha"] * deg2rad
            else:
                raise RuntimeError(f"Alpha must be in degrees or radians. Found {input_data['alpha'].unit} unit.")
        elif isinstance(input_data['alpha'], numbers.Real) or isinstance(input_data['alpha'], np.ndarray):
            if output_type == u.degree:
                output_dict['alpha'] = input_data['alpha'] * u.degree
            else:
                output_dict["alpha"] = np.deg2rad(input_data['alpha']) * u.radian
        else:
            raise RuntimeError("Alpha must be numeric with or without a unit.")
    return output_dict, found_radians


def check_alpha_requirement(path: List[str]):
    requires_alpha = False
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i + 1]
        requires_alpha = transform_graph.get_edge_data(step_start, step_end)['requires_alpha'] or requires_alpha
    return requires_alpha


def determine_input_kind(input_data: Dict[str, np.ndarray]) -> str:
    input_keys = list(input_data.keys())
    for valid_kind, param_list in VALID_KINDS.items():
        for param_option in param_list:
            if set(input_keys) == set(param_option):
                return valid_kind
    numeric_key_count = sum([isinstance(key, Quantity) for key in input_data])
    if numeric_key_count == 3:
        return "MZP"
    raise ValueError("Invalid Keys")  # TODO: improve this message


def get_transform_path(input_kind: str, output_kind: str) -> List[str]:
    try:
        path = nx.shortest_path(transform_graph, input_kind, output_kind)
    except nx.exception.NetworkXNoPath:
        raise RuntimeError(f"Not possible to convert {input_kind} to {output_kind}")  # TODO: make this a custom error
    return path


def get_transform_equation(path: List[str]) -> Callable:
    current_function = identity
    for i, step_start in enumerate(path[:-1]):
        step_end = path[i+1]
        current_function = _compose2(transform_graph.get_edge_data(step_start, step_end)['func'],
                                     current_function)
    return current_function


def _determine_image_shape(input_dict: Dict[str, np.ndarray]) -> Tuple[int, int]:
    keys = list(input_dict.keys())
    shape = input_dict[keys[0]].shape
    return shape


def add_alpha(input_data: Dict[str, np.ndarray], alpha_choice) -> Dict[str, np.ndarray]:
    # test if alpha exists. if not check if alpha keyword added. if not create default alpha with warning.

    img_shape = _determine_image_shape(input_data)
    if len(img_shape) == 2:  # it's an image and not just an array
        alpha_choice = "radial" if alpha_choice is None else alpha_choice
        if "alpha" not in input_data:
            if alpha_choice not in ALPHA_FUNCTIONS:
                raise ValueError(f"Requested a {alpha_choice} alpha type. "
                                 f"This is not valid. Must be in {ALPHA_FUNCTIONS.keys()}")
            input_data['alpha'] = ALPHA_FUNCTIONS[alpha_choice](img_shape)
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


def _compose2(f, g):
    return lambda *a, **kw: f(g(*a, **kw))


def identity(x):
    return x
