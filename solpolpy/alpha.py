"""Functions related to constructing an alpha array for transformation."""

import astropy.units as u
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord


def radial_north(shape):
    """An alpha array oriented north.

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
    - creates radial polarization map
    - angles increase in counterclockwise direction

    """
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.rot90(np.flipud(np.arctan2(yy, xx)), k=1))*u.radian


def radial_from_wcs(wcs, shape):
    """Construct a WCS-aware alpha angle field.

    This computes the per-pixel alpha angle from helioprojective coordinates,
    with solar north = 0 and counterclockwise-positive rotation. For
    partial-frame images, this samples the appropriate subset of the larger
    solar-centered alpha field instead of assuming the Sun is centered in the
    image array.
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
