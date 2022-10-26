from typing import Callable, Dict, List, Tuple
import functools
from collections import Counter
import numbers
from xml.dom.minidom import parseString

from astropy.io import fits
import astropy.units as u
from astropy.units.quantity import Quantity
import numpy as np
import networkx as nx
import warnings

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
        input_data = convert_image_list_to_dict(input_data, alpha=alpha)


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


def sanitize_data_dict(data_dict, output_angle_unit):
    if output_angle_unit not in [u.radian, u.degree]:
        raise RuntimeError(f"Keys must be converted to degrees or radians. Found output_angle_unit={output_angle_unit}")

    # Set up necessary conversions
    deg2rad = (np.pi * u.radian) / (180 * u.degree)
    rad2deg = (180 * u.degree) / (np.pi * u.radian)

    # Sanitize all the keys
    found_radians_key = False
    output_dict = dict()
    for key, value in data_dict.items():
        if isinstance(key, numbers.Real):
            if output_angle_unit == u.degree:
                output_dict[key * u.degree] = value
            else:
                output_dict[np.deg2rad(key) * u.radian] = value
        elif isinstance(key, Quantity):
            if key.unit not in [u.radian, u.degree]:
                raise RuntimeError(f"Unsupported key type of {key.unit} found. Must be radians or degrees.")
            if key.unit == u.radian:
                found_radians_key = True
                if output_angle_unit == u.degree:
                    output_dict[key*rad2deg] = value
                else:
                    output_dict[key] = value
            elif key.unit == u.degree:
                if output_angle_unit == u.radian:
                    output_dict[key*deg2rad] = value
                else:
                    output_dict[key] = value
        else:  # case for all the string keys
            output_dict[key] = value

    # Sanitize alpha map
    if "alpha" in data_dict:
        if isinstance(data_dict['alpha'], Quantity):
            if data_dict['alpha'].unit == u.radian:
                found_radians_key = True
                if output_angle_unit == u.degree:
                    output_dict["alpha"] = data_dict['alpha'] * rad2deg
                else:
                    output_dict["alpha"] = data_dict["alpha"]
            elif data_dict['alpha'].unit == u.degree:
                if output_angle_unit == u.degree:
                    output_dict["alpha"] = data_dict["alpha"]
                else:
                    output_dict["alpha"] = data_dict["alpha"] * deg2rad
            else:
                raise RuntimeError(f"Alpha must be in degrees or radians. Found {data_dict['alpha'].unit} unit.")
        elif isinstance(data_dict['alpha'], numbers.Real) or \
                (isinstance(data_dict['alpha'], np.ndarray) and np.issubdtype(data_dict['alpha'].dtype, np.number)):
            if output_angle_unit == u.degree:
                output_dict['alpha'] = data_dict['alpha'] * u.degree
            else:
                output_dict["alpha"] = np.deg2rad(data_dict['alpha']) * u.radian
        else:
            raise RuntimeError("Alpha must be numeric with or without a unit.")
    return output_dict, found_radians_key


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
    data_out={}

    for xlist_item in input_data:
        with fits.open(xlist_item) as hdul:
            if hdul[0].header['POLAR']== 'Clear':
                key_value='Clear'
            elif isinstance(hdul[0].header['POLAR'], numbers.Real):
                key_value=hdul[0].header['POLAR']*u.degree
            else:
                raise Exception("Didn't recognise the POLAR keyword")
            data_out[key_value]=hdul[0].data

    return data_out


def _convert_LASCO_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    # TODO: check that by converting the polar angle to 0->360 degree form doesn't require a modification to the dataframe (possible direction)
    # TODO: verify it's correct to add 360 degrees to -ve angles
    # TODO: verify if you reverese the angle (-60-> 120 & 60->240) the brightness is -ve.
    data_out={}


    for xlist_item in input_data:
        with fits.open(xlist_item) as hdul:
            if hdul[0].header['POLAR'] =='+60 Deg':
                key_value=60*u.degree 
            elif hdul[0].header['POLAR'] =='0 Deg':
                key_value=0*u.degree
            elif hdul[0].header['POLAR'] =='-60 Deg':
                key_value=-60*u.degree
                #key_value=120*u.degree
                #key_value=30*u.degree
                #key_value=300*u.degree
            elif hdul[0].header['POLAR'] =='Clear':
                key_value='Clear'
            else:
                raise Exception("Didn't recognise the POLAR keyword")

            data_out[key_value]=hdul[0].data

    return data_out

def convert_image_list_to_dict(input_data: List[str], alpha=None) -> Dict[str, np.ndarray]:
    # create output dictionary
    data_out={}
    
    # create list of FITS
    fits_type=[]

    # get length of list to determine how many files to process.
    list_len=len(input_data)
    assert list_len >= 2, 'requires at least 2 FITS files'

    for xlist_item in input_data:
        with fits.open(xlist_item) as hdul:
            fits_type.append(hdul[0].header['DETECTOR'])

    if len(set(fits_type)) != 1:
        raise Exception("Input FITS are of different types")

    if fits_type[0] == 'COR1' or fits_type[0] == 'COR2':
        data_out = _convert_STEREO_list_to_dict(input_data)
        alpha='radial90'
        # TODO: change alpha only if None, and also make a warning
    elif fits_type[0] == 'C2' or fits_type[0] == 'C3':
        data_out = _convert_LASCO_list_to_dict(input_data)
        alpha='radial'
        # TODO: change alpha only if None, and also make a warning
    else:
        raise Exception("the input FITS type is not supported. Use dictionary input.")

    # check if all polarized data entries, or if a B, pB pair
    if list_len==2:

        assert "Clear" in data_out, "expected a Clear (total brightness) data frame in the input polarization data pair"

        # checking that there is a polarized data frame, and what the value of the angle is.
        angle=None
        for key in data_out.keys():
            if isinstance(key, Quantity): 
                angle=key
        
        assert angle is not None, "the input polarization angle key is None, expected a numeric value"

        # to make BpB we need to add an alpha
        data_out = add_alpha(data_out, alpha)

        # resolve pB in terms of the radiance through a single arbitrary polarizer
        data_out = {'pB': pB_from_single_angle(data_out['Clear'], data_out[angle], angle, data_out['alpha']),
                     'B': data_out['Clear'], 
                     'alpha': data_out['alpha']}

    # check angles of keys - create warning if not appropriate angles.
    elif list_len>2:

        # check that the key angles are complimentary
        key_total=0*u.degree
        for key in data_out.keys():
            if key=='Clear':
                raise Exception("More than 3 files received. At least one input file was not polarized. Include a B, pB pair, or at least 3 complimentary polarized data frames." )
                
            key_total=key_total+key

        if key_total != 360*u.degree:
            warnings.warn("Input angles are not complimentary (total=360), processing but may not be accruate", Warning )

    return data_out


def _compose2(f, g):
    return lambda *a, **kw: f(g(*a, **kw))


def identity(x):
    return x


def pB_from_single_angle(B, B_theta, theta, alpha, tol=1e-6):
    """
    Converts unpolarized brightness,`B`, Radiance through a polarizer at angle theta,`B_theta`,
    Polarizer angle,`theta`, and Solar position angle of an image point, `alpha` into Coronal
    polarized brightness, `pB`.
    
    This function takes in four vars of `B`, `B_theta`, `theta`,and `alpha`.

    Parameters
    ----------
    B : np.ndarray
      'Clear' or total brightness data frame  
    B_theta : np.ndarray
      polarized data at angle theta
    theta : Quantity (astropy)
      angle of polarized data frame
    alpha : Quantity (astropy)
      alpha array to accompany STEREO or LASCO dataframes
    tol : float
      tolerence at which the denominator is converted to a nan (see Notes)

    Returns
    -------
    polarized brightness in terms of the radiance through a single arbitrary polarizer

    Notes
    ------
    see Equation 5 in Deforest et al. 2022.
    if this equation is used for single values, alpha will need to be a Quantity

    Equation (5) is problematic, because the denominator is small when theta-alpha 
    is near pm pi/4. Therefore, a nan is inserted when < tol

    Due to the cosine in the denominator, with an alpha varying from 0 to 360 degrees pB should vary from positive to
    negative twice crossing zero four times (pB will go to infinity at this point).
    """
    
    pB_denominator = np.cos( 2*(theta - alpha) )
    pB_denominator[np.abs(pB_denominator) < tol] = np.nan

    return ( B-( 2*B_theta ) ) / pB_denominator

