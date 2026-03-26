"""pytest test suite for the polarizers module of solpolpy
"""

import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from astropy.io import fits
from ndcube import NDCollection, NDCube

from solpolpy.alpha import radial_north
from solpolpy.core import (
    _determine_image_shape,
    add_alpha,
    determine_input_kind,
    get_transform_equation,
    get_transform_path,
    resolve,
)
from solpolpy.errors import UnsupportedTransformationError
from solpolpy.transforms import System
from solpolpy.util import solnorth_from_wcs
from tests.fixtures import *

# Solar WCS
wcs_sol = astropy.wcs.WCS(naxis=2)
wcs_sol.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
wcs_sol.wcs.cunit = "deg", "deg"
wcs_sol.wcs.cdelt = 0.5, 0.4
wcs_sol.wcs.crpix = 2, 2
wcs_sol.wcs.crval = 0.5, 1
wcs_sol.wcs.cname = "HPC lon", "HPC lat"

# Celestial WCS
wcs_cel = astropy.wcs.WCS(naxis=2)
wcs_cel.wcs.ctype = "RA---TAN", "DEC--TAN"
wcs_cel.wcs.cunit = "deg", "deg"
wcs_cel.wcs.cdelt = -0.5, 0.4
wcs_cel.wcs.crpix = 2.0, 2.0
wcs_cel.wcs.crval = 10.0, 20.0
wcs_cel.wcs.cname = "RA", "DEC"

hdr = fits.Header()
hdr.update(wcs_sol.to_header())
hdr.update(wcs_cel.to_header(key="A"))

wcs = astropy.wcs.WCS(hdr)


def test_determine_image_shape():
    data_out = [("B", NDCube(np.zeros((20, 20)), wcs=wcs, meta={"POLAR": "B"}))]
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert _determine_image_shape(collection) == (20, 20)


def test_add_alpha(bpb_ones_no_alpha):
    data_out = [("B", NDCube(np.zeros((10, 10)), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.zeros((10, 10)), wcs=wcs, meta={"POLAR": "pB"}))]
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert "alpha" not in collection
    out = add_alpha(collection)
    assert "alpha" in out


@pytest.mark.parametrize(
    "fixture_name, kind",
    [("npol_ones", System.npol),
     ("fourpol_ones", System.fourpol),
     ("mzpsolar_ones", System.mzpsolar),
     ("mzpsolar_ones_alpha", System.mzpsolar),
     ("bpb_ones", System.bpb),
     ("bpb_ones_no_alpha", System.bpb),
     ("btbr_ones", System.btbr),
     ("btbr_ones_no_alpha", System.btbr),
     ("stokes_ones", System.stokes),
     ("bp3_ones", System.bp3),
     ("bp3_ones_no_alpha", System.bp3),
     ("bthp_ones", System.bthp),
     ("mzpinstru_ones", System.mzpinstru)],
)
def test_determine_input_kind(fixture_name, kind, request):
    assert determine_input_kind(request.getfixturevalue(fixture_name)) == kind


def test_determine_input_kind_fail(example_fail):
    with pytest.raises(UnsupportedTransformationError):
        determine_input_kind(example_fail)


def test_check_all_paths_resolve(request):
    for source_system in System:
        for target_system in System:
            print(source_system, target_system)
            if source_system != target_system:  # don't need to convert if the input and output systems are the same
                try:
                    path = get_transform_path(source_system, target_system)
                except UnsupportedTransformationError:
                    pass  # skip any time the transform shouldn't be possible
                else:
                    print(source_system, target_system, path)
                    kwargs = {"out_angles": [0] * u.deg}
                    if getattr(get_transform_equation(path), "uses_in_angles", False):
                        kwargs["in_angles"] = [0, 60, 120] * u.deg
                    result = resolve(request.getfixturevalue(f"{source_system}_ones"),
                                     target_system,
                                     **kwargs)
                    assert isinstance(result, NDCollection)


def test_btbr_to_npol_missing_out_angles(btbr_ones):
    with pytest.raises(ValueError):
        resolve(btbr_ones, "npol")


def test_mzp_to_mzp_is_constant(mzp_ones_other_order):
    result = resolve(mzp_ones_other_order, "mzpsolar")
    assert result == mzp_ones_other_order


def test_mzp_to_npol_as_mzp_is_constant(mzpsolar_ones):
    result = resolve(mzpsolar_ones, "npol", out_angles=[-60, 0, 60] * u.degree)
    assert np.allclose(result['60.0 deg'].data, mzpsolar_ones["P"].data)
    assert np.allclose(result['-60.0 deg'].data, mzpsolar_ones["M"].data)
    assert np.allclose(result['0.0 deg'].data, mzpsolar_ones["Z"].data)


def test_add_alpha_uses_wcs_information():
    collection = NDCollection(
        [
            ("B", NDCube(np.zeros((6, 6)), wcs=wcs, meta={"POLAR": "B"})),
            ("pB", NDCube(np.zeros((6, 6)), wcs=wcs, meta={"POLAR": "pB"})),
        ],
        meta={},
        aligned_axes="all",
    )

    out = add_alpha(collection)

    assert "alpha" in out
    assert out["alpha"].data.shape == (6, 6)
    assert not np.allclose(u.Quantity(out["alpha"].data, unit=u.radian).value, 0)


def test_add_alpha_is_north_referenced_for_north_up_wcs():
    north_up_wcs = astropy.wcs.WCS(naxis=2)
    north_up_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    north_up_wcs.wcs.cunit = "deg", "deg"
    north_up_wcs.wcs.cdelt = 0.2, -0.2
    north_up_wcs.wcs.crpix = 3, 3
    north_up_wcs.wcs.crval = 0, 0

    collection = NDCollection(
        [
            ("B", NDCube(np.zeros((5, 5)), wcs=north_up_wcs, meta={"POLAR": "B"})),
            ("pB", NDCube(np.zeros((5, 5)), wcs=north_up_wcs, meta={"POLAR": "pB"})),
        ],
        meta={},
        aligned_axes="all",
    )

    out = add_alpha(collection)
    alpha = u.Quantity(out["alpha"].data, unit=u.radian).to_value(u.radian)

    np.testing.assert_allclose(alpha[0, 2], 0.0, atol=5e-3)
    np.testing.assert_allclose(alpha[2, 0], np.pi / 2, atol=5e-3)
    np.testing.assert_allclose(alpha[2, 4], -np.pi / 2, atol=5e-3)


def test_add_alpha_uses_wcs_for_off_center_partial_frame():
    off_center_wcs = astropy.wcs.WCS(naxis=2)
    off_center_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    off_center_wcs.wcs.cunit = "deg", "deg"
    off_center_wcs.wcs.cdelt = 0.2, -0.2
    off_center_wcs.wcs.crpix = 3, 3
    off_center_wcs.wcs.crval = 1.0, 0.0

    collection = NDCollection(
        [
            ("B", NDCube(np.zeros((5, 5)), wcs=off_center_wcs, meta={"POLAR": "B"})),
            ("pB", NDCube(np.zeros((5, 5)), wcs=off_center_wcs, meta={"POLAR": "pB"})),
        ],
        meta={},
        aligned_axes="all",
    )

    out = add_alpha(collection)
    alpha = u.Quantity(out["alpha"].data, unit=u.radian).to_value(u.radian)

    np.testing.assert_allclose(alpha[2, 2], -np.pi / 2, atol=5e-3)
    assert not np.allclose(alpha, radial_north((5, 5)).to_value(u.radian))


def test_solnorth_from_wcs_uses_north_zero_ccw_positive():
    north_up_wcs = astropy.wcs.WCS(naxis=2)
    north_up_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    north_up_wcs.wcs.cunit = "deg", "deg"
    north_up_wcs.wcs.cdelt = 0.2, -0.2
    north_up_wcs.wcs.crpix = 3, 3
    north_up_wcs.wcs.crval = 0, 0

    north_up = solnorth_from_wcs(north_up_wcs, (5, 5)).to_value(u.deg)
    np.testing.assert_allclose(north_up[2, 2], 0.0, atol=5e-3)

    rotated_wcs = astropy.wcs.WCS(naxis=2)
    rotated_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    rotated_wcs.wcs.cunit = "deg", "deg"
    rotated_wcs.wcs.cdelt = 0.2, -0.2
    rotated_wcs.wcs.crpix = 3, 3
    rotated_wcs.wcs.crval = 0, 0
    rotated_wcs.wcs.pc = np.array([[0, 1], [-1, 0]])

    rotated = solnorth_from_wcs(rotated_wcs, (5, 5)).to_value(u.deg)
    np.testing.assert_allclose(rotated[2, 2], 90.0, atol=5e-3)


def test_solnorth_from_wcs_matches_latitude_gradient_direction():
    north_up_wcs = astropy.wcs.WCS(naxis=2)
    north_up_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    north_up_wcs.wcs.cunit = "deg", "deg"
    north_up_wcs.wcs.cdelt = 0.2, -0.2
    north_up_wcs.wcs.crpix = 3, 3
    north_up_wcs.wcs.crval = 0, 0

    shape = (5, 5)
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    angle = solnorth_from_wcs(north_up_wcs, shape).to_value(u.rad)

    coords = astropy.wcs.utils.pixel_to_skycoord(x, y, north_up_wcs)
    lat = coords.Ty.to_value(u.deg)
    dy_lat, dx_lat = np.gradient(lat)
    norm = np.hypot(dx_lat, dy_lat)
    norm[norm == 0] = np.nan

    vx_from_angle = np.sin(angle)
    vy_from_angle = np.cos(angle)

    np.testing.assert_allclose(vx_from_angle, dx_lat / norm, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(vy_from_angle, -dy_lat / norm, atol=1e-12, rtol=1e-12)
