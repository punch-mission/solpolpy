import astropy.units as u
import numpy as np
from astropy.wcs import WCS
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
        crota = np.arccos(wcs.wcs.cd[0, 0]/cdelt1)
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
        return np.arctan2(wcs.wcs.pc[1, 0]/delta_ratio, wcs.wcs.pc[0, 0]) * u.rad
    elif hasattr(wcs.wcs, 'cd'):
        new_wcs = convert_cd_matrix_to_pc_matrix(wcs)
        return extract_crota_from_wcs(new_wcs)
    else:
        return 0 * u.rad
