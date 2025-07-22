"""Instrument specific code."""

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from ndcube import NDCollection, NDCube

from solpolpy.errors import TooFewFilesError, UnsupportedInstrumentError


def get_data_angle(header):
    angle_hdr = header["POLAR"]
    if isinstance(angle_hdr, float):
        angle = angle_hdr * u.degree
    elif isinstance(angle_hdr, str):
        angle = float(angle_hdr.split("Deg")[0]) * u.degree
    else:
        try:
            angle = float(angle_hdr) * u.degree
        except ValueError:
            msg = "Polar angle in the header could not be read for this instrument."
            raise UnsupportedInstrumentError(msg)

    return angle


def load_data(path_list: list[str],
              mask: np.ndarray | None = None,
              use_instrument_mask: bool = False,
              hdu_index: int = 0) -> NDCollection:
    """Basic loading function. See `load_with_occulter_mask`.

    Parameters
    ----------
    path_list: List[str]
        list of paths to FITS files to be loaded

    mask: Optional[np.ndarray]
        An optional mask to couple with images. Overrides any mask computed when use_instrument_mask=True.

    use_instrument_mask: bool
        If true, loads an instrument mask for common instruments defined in `get_instrument_mask`.

    hdu_index: int
        The index of the HDU (zero-based) to load data from.
        For many spacecraft this should be 0. For PUNCH it should be set to 1.

    Returns
    -------
    NDCollection
        The data are loaded as NDCollection object with WCS and header information available.
        The keys are labeled as 'angle_1', 'angle_2, 'angle_3', ...

    """
    # get length of list to determine how many files to process.
    if len(path_list) < 2:
        msg = "Requires at least 2 FITS files"
        raise TooFewFilesError(msg)

    data_out = []
    for _i, data_path in enumerate(path_list):
        with fits.open(data_path) as hdul:
            wcs = WCS(hdul[hdu_index].header)

            if use_instrument_mask and mask is None:
                mask = get_instrument_mask(hdul[hdu_index].header)

            if mask is None:  # make a mask of False if none is provided
                mask = np.zeros(hdul[hdu_index].data.shape, dtype=bool)

            angle = get_data_angle(hdul[hdu_index].header)

            data_out.append((str(angle),
                             NDCube(hdul[hdu_index].data,
                                    mask=mask,
                                    wcs=wcs,
                                    meta=dict(hdul[hdu_index].header))))

    return NDCollection(data_out, meta={}, aligned_axes="all")


def construct_mask(inner_radius: float,
                   outer_radius: float,
                   center_x: int,
                   center_y: int,
                   shape: tuple[int, int]) -> np.ndarray:
    """Constructs a mask where False indicates the pixel is valid and True masks the pixel as invalid.

    Pixels with radial distance between `inner_radius` and `outer_radius` are valid. Every other pixel is invalid.
    This is a standard coronograph mask.

    Parameters
    ----------
    inner_radius : float
        inner radius of the mask. pixels closer to the center than this are masked as invalid
    outer_radius : float
        outer radius of the mask. pixels farther from the center than this are masked as invalid
    center_x : int
        center pixel x coordinate
    center_y : int
        center pixel y coordinate
    shape : Tuple[int, int]
        the image dimensions

    Returns
    -------
    np.ndarray
        a coronograph mask that marks invalid pixels
        (those with radius less than `inner_radius` or greater than `outer_radius`) as True

    """
    xx, yy = np.ogrid[0:shape[0], 0:shape[1]]
    mask = np.zeros(shape, dtype=bool)
    mask[(xx - center_x) ** 2 + (yy - center_y) ** 2 < inner_radius] = True
    mask[(xx - center_x) ** 2 + (yy - center_y) ** 2 > outer_radius] = True
    return mask


def get_instrument_mask(header: fits.Header) -> np.ndarray:
    """Gets a coronograph mask for common instruments.

    Supports common solar instruments: LASCO, COSMO K, SECCHI

    Parameters
    ----------
    header : astropy.io.fits.Header
        a FITS header for the file in question

    Returns
    -------
    np.ndarray
        a pixel mask where valid pixels are marked False, masked out invalid pixels are flagged as True

    """
    x_shape = header["NAXIS1"]
    y_shape = header["NAXIS2"]

    center_x = header["CRPIX1"]
    center_y = header["CRPIX2"]

    R_sun = 960. / header["CDELT1"]  # TODO: clarify where 960. comes from

    radius = 0
    if header["INSTRUME"] == "LASCO":
        if "C2" in  header["DETECTOR"]:
            radius = 2.5
        elif "C3" in header["DETECTOR"]:
            radius = 4.0
        else:
            msg = "The data in this file is not formatted according to a known instrument."
            raise UnsupportedInstrumentError(msg)
    elif header["INSTRUME"] == "COSMO K-Coronagraph":
        radius = 1.15
    elif header["INSTRUME"] == "SECCHI":
        if "COR1" in header["DETECTOR"]:
            radius = 1.57
        elif "COR2" in header["DETECTOR"]:
            radius = 3.0
        else:
            msg = "The data in this file is not formatted according to a known instrument."
            raise UnsupportedInstrumentError(msg)
    else:
        msg = "The data in this file is not formatted according to a known instrument."
        raise UnsupportedInstrumentError(msg)

    y_array = [(j_step - center_y) / R_sun for j_step in range(y_shape)]
    outer_edge = np.max(y_array) * R_sun

    return construct_mask((radius * R_sun) ** 2, outer_edge ** 2, center_x, center_y, (x_shape, y_shape))
