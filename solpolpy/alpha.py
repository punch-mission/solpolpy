import numpy as np
import astropy.units as u


def radial(shape):
    '''
    assumes solar north is up
    creates radial polarization map
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.rot90(np.arctan2(yy, xx), k=1))*u.radian


ALPHA_FUNCTIONS = {'radial': radial}
