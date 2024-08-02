"""Functions related to constructing an alpha array for transformation."""

import astropy.units as u
import numpy as np


def radial_north(shape):
    """An alpha array oriented west.

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
    return np.rot90(np.fliplr(np.arctan2(yy, xx)+np.pi), k=1)*u.radian


ALPHA_FUNCTIONS = {"radial_north": radial_north,
                   "zeros": np.zeros}
