"""Instrument specific code"""
import warnings
import typing as t

from ndcube import NDCube, NDCollection
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


def load_data(path_list: t.List[str]) -> NDCollection:
    """
    Parameters
    ----------
    path_list: List[str] 
        list of paths to be loaded

    Returns
    -------
    NDCollection
        The data are loaded as NDCollection object with WCS and header information available.
        The keys are labeled as 'angle_1', 'angle_2, 'angle_3', ...
    """
    # get length of list to determine how many files to process.
    list_len = len(path_list)
    assert list_len >= 2, 'requires at least 2 FITS files'
    
    data_out = []
    for i, data_path in enumerate(path_list):
        with fits.open(data_path) as hdul:
            wcs = WCS(hdul[0].header)
            data_out.append(("angle_" + str(i+1),
                             NDCube(hdul[0].data, 
                                    wcs=wcs, 
                                    meta=hdul[0].header)))
            
    return NDCollection(data_out, meta={}, aligned_axes="all")


def load_with_occulter_mask(file_path: str,
                            mask_radii: t.Optional[float] = None,
                            center_pix_x: t.Optional[int] = None,
                            center_pix_y: t.Optional[int] = None) -> NDCube:
    """Creates a simple masked dataset for coronagraph data

    Parameters
    -----
    file_path : str
        path to a FITS file to be masked

    mask_radii : Optional[float]
        specifies the size of the mask in Solar Radii.
            If provided overwrites default mask sizes dependent on instruments.

    center_pix_x : Optional[int]
        if defined specifies the x position center of the masked region. If not provided, uses image center.

    center_pix_y : Optional[int]
        if defined specifies the y position center of the masked region. If not provided, uses image center.

    Returns
    ------
    NDcube
        An ndcube of output data, wcs and mask
    """
    file_of_interest = fits.open(file_path)
    file_data = file_of_interest[0].data
    file_meta = file_of_interest[0].header

    if file_meta['CRPIX1']:
        center_x = file_meta['CRPIX1']
    else:
        center_x = int(file_data.shape[0] / 2)
        warnings.warn("No CRPIX1 found in FITS header, setting CRPIX1 to center pixel")

    if file_meta['CRPIX2']:
        center_y = file_meta['CRPIX2']
    else:
        center_y = int(file_data.shape[1] / 2)
        warnings.warn("No CRPIX2 found in FITS header, setting CRPIX2 to center pixel")

    if file_meta['CDELT1']:
        R_sun = 960. / file_meta['CDELT1']  # this needs to be fixed with explanation, always the same?
    else:
        R_sun = 50
        warnings.warn("No CDELT1 found in FITS header, setting solar_radii to 50")

    if file_meta['INSTRUME'] == 'LASCO':
        detector_name = file_meta['DETECTOR']
        if 'C2' in detector_name:
            mask_R = 2.5
        elif 'C3' in detector_name:
            mask_R = 4.0
    elif file_meta['INSTRUME'] == 'COSMO K-Coronagraph':
        mask_R = 1.15
    elif file_meta['INSTRUME'] == 'SECCHI':
        detector_name = file_meta['DETECTOR']
        if 'COR1' in detector_name:
            mask_R = 1.57
        elif 'COR2' in detector_name:
            mask_R = 3.0
    else:
        warnings.warn("No INSTRUME found in FITS header, setting mask_radii to 1")
        mask_R = 1.0

    # user defined inputs to overwrite pre-defined inputs
    if mask_radii is not None:
        mask_R = mask_radii

    if center_pix_x is not None:
        center_x = center_pix_x

    if center_pix_y is not None:
        center_y = center_pix_y

    # obtain array dimensions
    x_shape, y_shape = file_data.shape

    # Calibrate to solar radius
    x_array = np.empty(x_shape)
    y_array = np.empty(y_shape)

    for i_step in range(x_shape):
        x_array[i_step] = (i_step - center_x) / R_sun

    for j_step in range(y_shape):
        y_array[j_step] = (j_step - center_y) / R_sun

    rr = np.max(y_array) * R_sun
    xx, yy = np.ogrid[0:x_shape, 0:y_shape]
    output_mask = np.ones(file_data.shape)  # the final mask array
    mask = (xx - center_x) ** 2 + (yy - center_y) ** 2 < (mask_R * R_sun) ** 2
    output_mask[mask] = 0
    mask = (xx - center_x) ** 2 + (yy - center_y) ** 2 > rr ** 2
    output_mask[mask] = 0

    return NDCube(data=file_data, wcs=WCS(file_meta), meta=file_meta, mask=output_mask)
