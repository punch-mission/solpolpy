import copy as copy

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.wcs import WCS, DistortionLookupTable
from astropy.wcs.utils import proj_plane_pixel_scales
from ndcube import NDCollection


def combine_all_collection_masks(collection: NDCollection) -> np.ndarray | None:
    """Combine all the masks in a given collection."""
    return combine_masks(*[cube.mask for key, cube in collection.items() if key != "alpha"])


def combine_masks(*args) -> np.ndarray | None:
    """Combine masks.

    If any of the masks are None, the result is None.
    Otherwise, when combining any value that is masked in any of the input args, gets masked, i.e. it does a logical or.
    """
    if any(arg is None for arg in args):
        return None
    else:
        return np.logical_or.reduce(args)


def calculate_pc_matrix(crota: float, cdelt: (float, float)) -> np.ndarray:
    """
    Calculate a PC matrix given CROTA and CDELT.

    Parameters
    ----------
    crota : float
        rotation angle from the WCS
    cdelt : float
        pixel size from the WCS

    Returns
    -------
    np.ndarray
        PC matrix

    """
    return np.array(
        [
            [np.cos(crota), -np.sin(crota) * (cdelt[0] / cdelt[1])],
            [np.sin(crota) * (cdelt[1] / cdelt[0]), np.cos(crota)],
        ],
    )


def convert_cd_matrix_to_pc_matrix(wcs):
    if hasattr(wcs.wcs, 'cd'):
        cdelt1, cdelt2 = proj_plane_pixel_scales(wcs)
        crota = np.arccos(wcs.wcs.cd[0, 0] / cdelt1)
        new_wcs = WCS(naxis=2)
        new_wcs.wcs.ctype = wcs.wcs.ctype
        new_wcs.wcs.crval = wcs.wcs.crval
        new_wcs.wcs.crpix = wcs.wcs.crpix
        new_wcs.wcs.pc = calculate_pc_matrix(crota, (cdelt1, cdelt2))
        new_wcs.wcs.cdelt = (-cdelt1, cdelt2)
        new_wcs.wcs.cunit = 'deg', 'deg'
        return new_wcs
    else:  # noqa RET505
        return wcs


def extract_crota_from_wcs(wcs: WCS) -> u.deg:
    """Extract CROTA from a WCS."""
    if hasattr(wcs.wcs, 'pc'):
        delta_ratio = wcs.wcs.cdelt[1] / wcs.wcs.cdelt[0]
        return np.arctan2(wcs.wcs.pc[1, 0] / delta_ratio, wcs.wcs.pc[0, 0]) * u.rad
    elif hasattr(wcs.wcs, 'cd'):
        new_wcs = convert_cd_matrix_to_pc_matrix(wcs)
        return extract_crota_from_wcs(new_wcs)
    else:
        return 0 * u.rad


def indexed_offset(index, distortion, image_shape):
    return distortion.get_offset(index % image_shape[0], index // image_shape[0])


def calculate_distortion(distortion, image_shape):
    vectorized_offset = np.vectorize(indexed_offset, excluded=["distortion", "image_shape"])
    distortion_array = vectorized_offset(np.arange(image_shape[0] * image_shape[1]),
                                         distortion=distortion,
                                         image_shape=image_shape)
    return distortion_array.reshape(image_shape)


def compute_distortion_shift(image_shape, wcs: WCS):
    """
    Calculate shift in pixels due to optical distortion

    Parameters
    ----------
    image_shape : Tuple(int, int)
        shape of input image
    wcs : WCS
        WCS from input object

    Returns
    -------
    tuple (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray)
         Tuple containing new x-coordinates, new y-coordinates, valid mask,
         original i-coordinates, and original j-coordinates.

    """
    shift_x = calculate_distortion(wcs.cpdis1, image_shape)
    shift_y = calculate_distortion(wcs.cpdis2, image_shape)
    # Calculate new coordinates for pixel shifts
    i_coords, j_coords = np.meshgrid(np.arange(image_shape[0]), np.arange(image_shape[1]), indexing='ij')
    new_x = np.round(j_coords + shift_x).astype(int)
    new_y = np.round(i_coords + shift_y).astype(int)
    valid_mask = (0 <= new_x) & (new_x < image_shape[1]) & (0 <= new_y) & (new_y < image_shape[0])

    return new_x, new_y, valid_mask, i_coords, j_coords


def apply_distortion_shift(input_image, new_x, new_y, valid_mask, i_coords, j_coords):
    """
    Apply shift in pixels due to optical distortion

    Parameters
    ----------
    input_image : np.ndarray
        input image on which shift to be applied
    new_x : np.ndarray
        The precomputed new x-coordinates after distortion.
    new_y : np.ndarray
        The precomputed new y-coordinates after distortion.
    valid_mask : np.ndarray
        Boolean mask indicating valid shifts within bounds.
    i_coords : np.ndarray
        Original i-coordinates of pixels frpm input_image.
    j_coords : np.ndarray
        Original j-coordinates of pixels from input_image.
    Returns
    -------
    np.ndarray
        Image after applying the distortion shifts.
    """
    shifted_image = copy.copy(input_image)
    shifted_image[new_y[valid_mask], new_x[valid_mask]] = input_image[i_coords[valid_mask], j_coords[valid_mask]]
    return shifted_image


def make_empty_distortion_model(num_bins: int, image: np.ndarray) -> (DistortionLookupTable, DistortionLookupTable):
    """ Create an empty distortion table

    Parameters
    ----------
    num_bins : int
        number of histogram bins in the distortion model, i.e. the size of the distortion model is (num_bins, num_bins)
    image : np.ndarray
        image to create a distortion model for

    Returns
    -------
    (DistortionLookupTable, DistortionLookupTable)
        x and y distortion models
    """
    # make an initial empty distortion model
    r = np.linspace(0, image.shape[0], num_bins + 1)
    c = np.linspace(0, image.shape[1], num_bins + 1)
    r = (r[1:] + r[:-1]) / 2
    c = (c[1:] + c[:-1]) / 2

    err_px, err_py = r, c
    err_x = np.ones((num_bins, num_bins))
    err_y = np.zeros((num_bins, num_bins))

    cpdis1 = DistortionLookupTable(
        -err_x.astype(np.float32), (0, 0), (err_px[0], err_py[0]), ((err_px[1] - err_px[0]), (err_py[1] - err_py[0]))
    )
    cpdis2 = DistortionLookupTable(
        -err_y.astype(np.float32), (0, 0), (err_px[0], err_py[0]), ((err_px[1] - err_px[0]), (err_py[1] - err_py[0]))
    )
    return cpdis1, cpdis2


def collection_to_maps(collection):
    """
    Convert an NDCollection to a list of SunPy Map objects.

    Parameters:
    -----------
    ndcollection : NDCollection
        The NDCollection containing data, metadata and wcs.

    Returns:
    --------
    list of sunpy.map.Map
        A list of SunPy Map objects created from the NDCollection.
    """
    sunpy_maps = []

    for key in collection.keys():
        data = collection[key].data
        wcs = collection[key].wcs

        sunpy_maps.append(sunpy.map.Map(data, wcs))

    return sunpy_maps
