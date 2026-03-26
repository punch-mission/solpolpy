"""Functions related to constructing an alpha array for transformation."""

import astropy.units as u
import numpy as np
import sunpy.coordinates  # noqa: F401
from astropy.wcs.utils import pixel_to_skycoord


def radial_north(shape):
    """An alpha array referenced to north with counterclockwise-positive angles.

    Parameters
    ----------
    shape : tuple[int, int]
        how big the array should be

    Returns
    -------
    np.ndarray
        alpha array used in calculations

    Notes
    -----
    - assumes solar north is up
    - assumes polarizer 0 is along solar north axis
    - uses NumPy image indexing: row 0 is the top of the image, column 0 is the left
    - returns the radial axis angle measured from north = 0, increasing counterclockwise

    """
    nrows, ncols = shape
    center_row = (nrows - 1) / 2.0
    center_col = (ncols - 1) / 2.0

    row_indices, col_indices = np.indices(shape, dtype=float)
    dx = col_indices - center_col
    dy_up = center_row - row_indices

    # Angle from north with counterclockwise-positive rotation.
    return np.arctan2(-dx, dy_up) * u.radian


def radial_from_wcs(wcs, shape):
    """Construct an alpha array from solar coordinates in the WCS.

    This computes the radial direction from solar center for each pixel,
    measured from solar north = 0 with counterclockwise-positive rotation.
    For partial-frame images, this samples the relevant subset of the full
    solar-centered alpha field instead of assuming the Sun is at image center.
    """
    nrows, ncols = shape
    row_indices, col_indices = np.mgrid[0:nrows, 0:ncols]
    coords = pixel_to_skycoord(col_indices, row_indices, wcs)

    tx = coords.Tx.to_value(u.deg)
    ty = coords.Ty.to_value(u.deg)

    return np.arctan2(-tx, ty) * u.radian


ALPHA_FUNCTIONS = {"radial_north": radial_north,
                   "radial_from_wcs": radial_from_wcs,
                   "zeros": np.zeros}
