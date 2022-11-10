import numbers
from typing import List, Dict

import numpy as np
from astropy import units as u
from astropy.io import fits


def _convert_LASCO_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    # TODO: check that by converting the polar angle to 0->360 degree form doesn't require a modification to the dataframe (possible direction)
    # TODO: verify it's correct to add 360 degrees to -ve angles
    # TODO: verify if you reverese the angle (-60-> 120 & 60->240) the brightness is -ve.
    data_out = {}

    for xlist_item in input_data:
        with fits.open(xlist_item) as hdul:
            if hdul[0].header['POLAR'] == '+60 Deg':
                key_value = 60*u.degree
            elif hdul[0].header['POLAR'] == '0 Deg':
                key_value = 0*u.degree
            elif hdul[0].header['POLAR'] == '-60 Deg':
                key_value = -60*u.degree
            elif hdul[0].header['POLAR'] == 'Clear':
                key_value = 'Clear'
            else:
                raise Exception("Didn't recognise the POLAR keyword")

            data_out[key_value] = hdul[0].data

    return data_out


def _convert_STEREO_list_to_dict(input_data: List[str]) -> Dict[str, np.ndarray]:
    data_out = {}

    for xlist_item in input_data:
        with fits.open(xlist_item) as hdul:
            if hdul[0].header['POLAR'] == 'Clear':
                key_value = 'Clear'
            elif isinstance(hdul[0].header['POLAR'], numbers.Real):
                key_value = hdul[0].header['POLAR']*u.degree  # TODO: make sure it is -60, 0, 60 from now on
            else:
                raise Exception("Didn't recognise the POLAR keyword")
            data_out[key_value] = hdul[0].data
    return data_out
