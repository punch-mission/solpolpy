import numpy as np


def radial(shape):
    x_size, y_size = shape
    x = np.arange(-x_size // 2, x_size // 2)
    y = np.arange(-y_size // 2, y_size // 2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.rot90(np.arctan2(yy, xx), k=1))


ALPHA_FUNCTIONS = {'radial': radial}
