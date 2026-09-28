import astropy.units as u
import astropy.wcs
import numpy as np
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
    alpha_mask = np.ones((2, 2), dtype=bool)
    collection = NDCollection(
        [
            ("B", NDCube(np.ones((2, 2)), wcs=solar_wcs, mask=mask_b)),
            ("pB", NDCube(np.ones((2, 2)), wcs=solar_wcs, mask=mask_pb)),
            ("alpha", NDCube(np.zeros((2, 2)) * u.deg, wcs=solar_wcs, mask=alpha_mask)),
        ],
        aligned_axes="all",
    )

    combined = combine_all_collection_masks(collection)

    np.testing.assert_array_equal(combined, mask_b | mask_pb)


def test_calculate_pc_matrix_for_quarter_turn_and_unequal_pixel_scales():
    pc_matrix = calculate_pc_matrix(np.pi / 2, (2.0, 1.0))

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
    assert convert_cd_matrix_to_pc_matrix(solar_wcs) is solar_wcs


def test_extract_crota_from_pc_matrix():
    rotated_wcs = solar_wcs.deepcopy()
    rotated_wcs.wcs.pc = calculate_pc_matrix(np.deg2rad(30), rotated_wcs.wcs.cdelt)

    np.testing.assert_allclose(extract_crota_from_wcs(rotated_wcs), 30 * u.deg)


def test_empty_distortion_model_has_one_pixel_x_shift():
    image = np.zeros((3, 4))
    distortion_x, distortion_y = make_empty_distortion_model(2, image)

    assert indexed_offset(0, distortion_x, image.shape) == -1
    np.testing.assert_allclose(calculate_distortion(distortion_x, image.shape), -1)
    np.testing.assert_allclose(calculate_distortion(distortion_y, image.shape), 0)


def test_compute_and_apply_distortion_shift():
    image = np.arange(12).reshape(3, 4)
    distortion_wcs = astropy.wcs.WCS(naxis=2)
    distortion_wcs.cpdis1, distortion_wcs.cpdis2 = make_empty_distortion_model(2, image)

    shift_coordinates = compute_distortion_shift(image.shape, distortion_wcs)
    shifted_image = apply_distortion_shift(image, *shift_coordinates)

    # The model moves every valid source pixel one column to the left. The
    # unfilled right edge retains its original value.
    expected = np.array([[1, 2, 3, 3], [5, 6, 7, 7], [9, 10, 11, 11]])
    np.testing.assert_array_equal(shifted_image, expected)
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
    np.testing.assert_allclose(latitude[1, 1], 0, atol=1e-12)
    np.testing.assert_allclose(latitude[:, 1], [-0.2, 0, 0.2, 0.4, 0.6], atol=3e-5)


def test_solnorth_from_wcs_points_along_detector_y():
    latitude = compute_lats(solar_wcs, (5, 5))

    angle = solnorth_from_wcs(solar_wcs, (5, 5))
    angle_with_precomputed_lats = solnorth_from_wcs(solar_wcs, (5, 5), precomputed_lats=latitude)

    assert angle.unit == u.deg
    np.testing.assert_allclose(angle.to_value(u.deg), 0, atol=1e-9)
    np.testing.assert_allclose(angle_with_precomputed_lats, angle, atol=1e-9 * u.deg)


def test_wrap_pm_pi_preserves_units_and_wraps_upper_boundary():
    angle = np.array([-540, -180, -10, 180, 540]) * u.deg

    wrapped = wrap_pm_pi(angle)

    assert wrapped.unit == u.deg
    np.testing.assert_allclose(wrapped.value, [-180, -180, -10, -180, -180])
