import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from ndcube import NDCollection, NDCube
from sunpy.map import GenericMap

from solpolpy.util import (
    apply_distortion_shift,
    calculate_distortion,
    calculate_pc_matrix,
    collection_to_maps,
    combine_all_collection_masks,
    combine_masks,
    compute_distortion_shift,
    compute_lats,
    convert_cd_matrix_to_pc_matrix,
    extract_crota_from_wcs,
    indexed_offset,
    make_empty_distortion_model,
    solnorth_from_wcs,
    wrap_pm_pi,
)

# Solar WCS
solar_wcs = astropy.wcs.WCS(naxis=2)
solar_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
solar_wcs.wcs.cunit = "deg", "deg"
solar_wcs.wcs.cdelt = 0.2, 0.2
solar_wcs.wcs.crpix = 2, 2
solar_wcs.wcs.crval = 0, 0
solar_wcs.wcs.cname = "HPC lon", "HPC lat"


def test_combine_masks_uses_logical_or():
    mask_a = np.array([[False, True], [False, False]])
    mask_b = np.array([[False, False], [True, False]])

    combined = combine_masks(mask_a, mask_b)

    expected = np.array([[False, True], [True, False]])
    np.testing.assert_array_equal(combined, expected)


def test_combine_masks_returns_none_when_any_mask_is_none():
    mask = np.zeros((2, 2), dtype=bool)

    assert combine_masks(mask, None) is None


def test_combine_all_collection_masks_ignores_alpha_mask():
    mask_b = np.array([[False, True], [False, False]])
    mask_pb = np.array([[False, False], [True, False]])
    alpha_mask = np.array([[True, False], [False, False]])
    collection = NDCollection(
        [
            ("B", NDCube(np.ones((2, 2)), wcs=solar_wcs, mask=mask_b)),
            ("pB", NDCube(np.ones((2, 2)), wcs=solar_wcs, mask=mask_pb)),
            ("alpha", NDCube(np.zeros((2, 2)) * u.deg, wcs=solar_wcs, mask=alpha_mask)),
        ],
        aligned_axes="all",
    )

    combined = combine_all_collection_masks(collection)

    # Pixel (0, 0) is masked only in alpha, so it must remain unmasked if the
    # alpha mask is excluded from the logical OR.
    assert alpha_mask[0, 0]
    assert not mask_b[0, 0] and not mask_pb[0, 0]
    assert not combined[0, 0]
    np.testing.assert_array_equal(combined, mask_b | mask_pb)


def test_calculate_pc_matrix_for_quarter_turn_and_unequal_pixel_scales():
    pc_matrix = calculate_pc_matrix(np.pi / 2, (2.0, 1.0))

    # The PC matrix is [[cos(theta), -sin(theta) CDELT1/CDELT2],
    #                   [sin(theta) CDELT2/CDELT1, cos(theta)]].
    # For theta = pi/2 and CDELT = (2, 1), cos(theta) = 0,
    # sin(theta) = 1, and the two scale ratios are 2 and 1/2.
    expected = np.array([[0.0, -2.0], [0.5, 0.0]])
    np.testing.assert_allclose(pc_matrix, expected, atol=1e-15)


def test_convert_cd_matrix_to_pc_matrix_preserves_wcs_geometry():
    cd_wcs = astropy.wcs.WCS(naxis=2)
    cd_wcs.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
    cd_wcs.wcs.crval = 1, 2
    cd_wcs.wcs.crpix = 3, 4
    cd_wcs.wcs.cd = np.array([[0.0, -0.2], [0.1, 0.0]])

    pc_wcs = convert_cd_matrix_to_pc_matrix(cd_wcs)

    np.testing.assert_allclose(pc_wcs.wcs.crval, cd_wcs.wcs.crval)
    np.testing.assert_allclose(pc_wcs.wcs.crpix, cd_wcs.wcs.crpix)
    np.testing.assert_allclose(pc_wcs.wcs.cdelt, [-0.1, 0.2])
    np.testing.assert_allclose(extract_crota_from_wcs(pc_wcs), -90 * u.deg)


def test_convert_cd_matrix_to_pc_matrix_returns_pc_wcs_unchanged():
    pc_wcs = convert_cd_matrix_to_pc_matrix(solar_wcs)

    np.testing.assert_allclose(pc_wcs.wcs.crval, solar_wcs.wcs.crval)
    np.testing.assert_allclose(pc_wcs.wcs.crpix, solar_wcs.wcs.crpix)
    np.testing.assert_allclose(pc_wcs.wcs.cdelt, solar_wcs.wcs.cdelt)
    np.testing.assert_allclose(pc_wcs.wcs.pc, solar_wcs.wcs.pc)
    assert list(pc_wcs.wcs.ctype) == list(solar_wcs.wcs.ctype)
    assert list(pc_wcs.wcs.cunit) == list(solar_wcs.wcs.cunit)
    assert list(pc_wcs.wcs.cname) == list(solar_wcs.wcs.cname)


def test_extract_crota_from_pc_matrix():
    rotated_wcs = solar_wcs.deepcopy()
    rotated_wcs.wcs.pc = calculate_pc_matrix(np.deg2rad(30), rotated_wcs.wcs.cdelt)

    np.testing.assert_allclose(extract_crota_from_wcs(rotated_wcs), 30 * u.deg)


def test_empty_distortion_model_has_zero_offsets():
    image = np.zeros((3, 4))
    distortion_x, distortion_y = make_empty_distortion_model(2, image)

    assert indexed_offset(0, distortion_x, image.shape) == 0
    np.testing.assert_allclose(calculate_distortion(distortion_x, image.shape), 0)
    np.testing.assert_allclose(calculate_distortion(distortion_y, image.shape), 0)


def test_compute_and_apply_distortion_shift():
    image = np.arange(12).reshape(3, 4)
    distortion_wcs = astropy.wcs.WCS(naxis=2)
    distortion_wcs.cpdis1, distortion_wcs.cpdis2 = make_empty_distortion_model(2, image)

    shift_coordinates = compute_distortion_shift(image.shape, distortion_wcs)
    shifted_image = apply_distortion_shift(image, *shift_coordinates)

    # A zero-distortion model leaves every pixel at its original coordinate.
    np.testing.assert_array_equal(shifted_image, image)
    np.testing.assert_array_equal(image, np.arange(12).reshape(3, 4))


def test_collection_to_maps_preserves_data():
    collection = NDCollection(
        [
            ("M", NDCube(np.zeros((2, 2)), wcs=solar_wcs)),
            ("P", NDCube(np.ones((2, 2)), wcs=solar_wcs)),
        ],
        aligned_axes="all",
    )

    maps = collection_to_maps(collection)

    assert len(maps) == 2
    assert all(isinstance(solar_map, GenericMap) for solar_map in maps)
    np.testing.assert_array_equal(maps[0].data, collection["M"].data)
    np.testing.assert_array_equal(maps[1].data, collection["P"].data)


def test_compute_lats_matches_simple_solar_wcs():
    latitude = compute_lats(solar_wcs, (5, 5))

    assert latitude.shape == (5, 5)
    assert latitude.unit == u.deg
    np.testing.assert_allclose(latitude[:, 1].to_value(u.deg), [-0.2, 0, 0.2, 0.4, 0.6], atol=3e-5)


def test_solnorth_from_wcs_points_along_detector_y():
    latitude = compute_lats(solar_wcs, (5, 5))

    angle = solnorth_from_wcs(solar_wcs, (5, 5))
    angle_with_precomputed_lats = solnorth_from_wcs(solar_wcs, (5, 5), precomputed_lats=latitude)
    angle_with_radian_lats = solnorth_from_wcs(solar_wcs, (5, 5), precomputed_lats=latitude.to(u.radian))

    assert angle.unit == u.deg
    np.testing.assert_allclose(angle.to_value(u.deg), 0, atol=1e-9)
    np.testing.assert_allclose(angle_with_precomputed_lats, angle, atol=1e-9 * u.deg)
    np.testing.assert_allclose(angle_with_radian_lats, angle, atol=1e-9 * u.deg)


def test_solnorth_from_wcs_rejects_unitless_precomputed_lats():
    latitude = compute_lats(solar_wcs, (5, 5))

    with pytest.raises(TypeError, match="astropy Quantity"):
        solnorth_from_wcs(solar_wcs, (5, 5), precomputed_lats=latitude.value)


def test_wrap_pm_pi_preserves_units_and_wraps_upper_boundary():
    angle = np.array([-540, -180, -10, 180, 540]) * u.deg

    wrapped = wrap_pm_pi(angle)

    assert wrapped.unit == u.deg
    # The interval is [-180, 180), so +180 is excluded and uses the equivalent
    # -180 representation. Every odd multiple of 180 degrees wraps there too.
    np.testing.assert_allclose(wrapped.value, [-180, -180, -10, -180, -180])
