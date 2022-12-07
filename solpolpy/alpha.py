import numpy as np
import astropy.units as u


# TODO: make sure these have up/solar north as the reference of 0 degrees


def radial90(shape):
    '''
    assumes solar north is to the left - useful for STEREO
    creates radial polarization map
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.rot90(np.arctan2(yy, xx), k=1))*u.radian


def radial_north(shape):
    '''
    assumes solar north is up
    creates radial polarization map
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    # return np.fliplr(np.arctan2(yy, xx))*u.radian
    return np.rot90(np.fliplr(np.arctan2(yy, xx) + np.pi), k=1) * u.radian


def radial_west(shape):
    '''
    assumes solar north is up
    creates radial polarization map
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.arctan2(yy, xx) + np.pi)*u.radian
    # return np.rot90(np.fliplr(np.arctan2(yy, xx) + np.pi), k=1) * u.radian
#

def zeros(shape):
    return np.zeros(shape)


ALPHA_FUNCTIONS = {'radial_north': radial_north,
                   'radial_west': radial_west,
                   'radial90': radial90,
                   'zeros': zeros}
