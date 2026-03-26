import astropy.units as u
import astropy.wcs
import numpy as np
import pytest

from solpolpy.alpha import radial_from_wcs, radial_north


@pytest.mark.parametrize("shape", [(1000, 1000), (500, 500), (500, 100)])
def test_radial_north(shape):
    alpha = radial_north(shape)
    assert isinstance(alpha, np.ndarray)
    assert alpha.shape == shape
    assert np.isclose(alpha.min(), -np.pi*u.radian, atol=0.01)
    assert np.isclose(alpha.max(), np.pi*u.radian, atol=0.01)


def test_radial_north_matches_python_image_indexing_and_ccw_sign():
    alpha = radial_north((5, 5)).to_value(u.radian)

    np.testing.assert_allclose(alpha[0, 2], 0.0, atol=1e-12)
    np.testing.assert_allclose(alpha[2, 0], np.pi / 2, atol=1e-12)
    np.testing.assert_allclose(alpha[2, 4], -np.pi / 2, atol=1e-12)
    np.testing.assert_allclose(np.abs(alpha[4, 2]), np.pi, atol=1e-12)


def test_radial_from_wcs_matches_centered_north_up_geometry():
    wcs = astropy.wcs.WCS(naxis=2)
    wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    wcs.wcs.cunit = "deg", "deg"
    wcs.wcs.cdelt = 0.2, -0.2
    wcs.wcs.crpix = 3, 3
    wcs.wcs.crval = 0, 0

    alpha = radial_from_wcs(wcs, (5, 5)).to_value(u.radian)

    np.testing.assert_allclose(alpha[0, 2], 0.0, atol=5e-3)
    np.testing.assert_allclose(alpha[2, 0], np.pi / 2, atol=5e-3)
    np.testing.assert_allclose(alpha[2, 4], -np.pi / 2, atol=5e-3)


def test_radial_from_wcs_handles_off_center_partial_frame():
    wcs = astropy.wcs.WCS(naxis=2)
    wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    wcs.wcs.cunit = "deg", "deg"
    wcs.wcs.cdelt = 0.2, -0.2
    wcs.wcs.crpix = 3, 3
    wcs.wcs.crval = 1.0, 0.0

    alpha = radial_from_wcs(wcs, (5, 5)).to_value(u.radian)

    np.testing.assert_allclose(alpha[2, 2], -np.pi / 2, atol=5e-3)
