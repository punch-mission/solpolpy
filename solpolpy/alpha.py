import numpy as np
import astropy.units as u


def radial_west(shape):
    '''
    assumes solar north is up
    assumes polarizer 0 is along positive horizontal axis
    creates radial polarization map
    angles increase in counterclockwise direction
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.flipud(np.fliplr(np.arctan2(yy, xx) + np.pi))*u.radian


def radial_north(shape):
    '''
    assumes solar north is up
    assumes polarizer 0 is along solar north axis
    creates radial polarization map
    angles increase in counterclockwise direction
    '''
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.rot90(np.fliplr(np.arctan2(yy, xx)+np.pi), k=1)*u.radian
    

def zeros(shape):
    return np.zeros(shape)


ALPHA_FUNCTIONS = {'radial_north': radial_north,
                   'radial_west': radial_west,
                   'zeros': zeros}
