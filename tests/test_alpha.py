import astropy.units as u
import numpy as np
import pytest

from solpolpy.alpha import radial_north


@pytest.mark.parametrize("shape", [(1000, 1000), (500, 500), (500, 100)])
def test_radial_north(shape):
    alpha = radial_north(shape)
    assert isinstance(alpha, np.ndarray)
    assert alpha.shape == shape
    assert np.isclose(alpha.min(), 0.0*u.radian, atol=0.01)
    assert np.isclose(alpha.max(), 2*np.pi*u.radian, atol=0.01)
