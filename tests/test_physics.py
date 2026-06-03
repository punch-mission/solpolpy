from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from ndcube import NDCollection, NDCube

from solpolpy.errors import SolpolpyError
from solpolpy.physics import (
    MZP_ANGLES,
    as_angle,
    bp3_from_analyzer_brightness,
    bp3_to_analyzer_brightness,
    calculate_angle_difference,
    clone_meta,
    get_data_keys,
    get_template_cube,
    project_three_polarizer_brightness,
    solve_three_polarizer_brightness,
    stack_data,
    wrap_linear_polarization,
    wrap_pm_pi,
)
from tests.fixtures import wcs


def _collection():
    return NDCollection(
        [
            ("B", NDCube(np.ones((2, 2)), wcs=wcs, meta={"POLAR": "B"})),
            ("pB", NDCube(np.full((2, 2), 0.25), wcs=wcs, meta={"POLAR": "pB"})),
            ("alpha", NDCube(np.zeros((2, 2)) * u.rad, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )


def test_collection_helpers_ignore_alpha_and_preserve_order():
    collection = _collection()

    assert get_data_keys(collection) == ["B", "pB"]
    assert get_template_cube(collection) is collection["B"]
    assert get_template_cube(collection, preferred_key="pB") is collection["pB"]
    np.testing.assert_equal(stack_data(collection, ["pB", "B"])[0], collection["pB"].data)


def test_clone_meta_deep_copies_and_updates_metadata():
    cube = SimpleNamespace(meta={"POLAR": "B", "nested": {"value": 1}})

    meta = clone_meta(cube, POLAR="pB")

    assert meta["POLAR"] == "pB"
    meta["nested"]["value"] = 2
    assert cube.meta["nested"]["value"] == 1


def test_as_angle_applies_default_unit_to_unitless_values():
    actual = as_angle([0, 90], default_unit=u.degree)

    assert actual.unit == u.degree
    np.testing.assert_allclose(actual.value, [0, 90])


def test_quantity_input_rejects_unitless_wrapped_angles():
    with pytest.raises(TypeError):
        wrap_pm_pi(np.pi)


def test_angle_wrapping_preserves_input_units():
    wrapped = wrap_pm_pi(190 * u.degree)
    linear = wrap_linear_polarization(100 * u.degree)

    assert wrapped.unit == u.degree
    assert linear.unit == u.degree
    np.testing.assert_allclose(wrapped.value, -170)
    np.testing.assert_allclose(linear.value, -80)


def test_calculate_angle_difference_wraps_quantity_inputs():
    actual = calculate_angle_difference(10 * u.degree, 350 * u.degree)

    assert actual.unit == u.degree
    np.testing.assert_allclose(actual.value, 20)


def test_bp3_analyzer_brightness_roundtrip_recovers_inputs():
    B = np.array([[10.0, 12.0], [8.0, 9.0]])
    pB = np.array([[1.0, 1.5], [0.5, 0.75]])
    pBp = np.array([[0.2, -0.3], [0.1, -0.2]])
    alpha = np.full((2, 2), 15.0) * u.degree

    brightness = bp3_to_analyzer_brightness(B, pB, pBp, alpha, MZP_ANGLES)
    recovered = bp3_from_analyzer_brightness(brightness, MZP_ANGLES, alpha)

    np.testing.assert_allclose(recovered[0], B, atol=1e-12)
    np.testing.assert_allclose(recovered[1], pB, atol=1e-12)
    np.testing.assert_allclose(recovered[2], pBp, atol=1e-12)


def test_solve_three_polarizer_brightness_supports_pixelwise_reference_angle():
    observed_angles = np.array([-60.0, 0.0, 60.0]) * u.degree
    solved_angles = np.array([0.0, 60.0, 120.0]) * u.degree
    reference_angle = np.array([[0.0, 5.0], [10.0, 15.0]]) * u.degree
    source = np.stack(
        [
            np.full((2, 2), 3.0),
            np.full((2, 2), 2.0),
            np.full((2, 2), 1.0),
        ],
        axis=0,
    )

    projected = project_three_polarizer_brightness(source, solved_angles, observed_angles, reference_angle)
    recovered = solve_three_polarizer_brightness(projected, observed_angles, solved_angles, reference_angle)

    np.testing.assert_allclose(recovered, source, atol=1e-12)


def test_solve_three_polarizer_brightness_raises_for_degenerate_angles():
    brightness = np.ones((3, 2, 2))
    duplicate_angles = np.array([0.0, 0.0, 60.0]) * u.degree

    with pytest.raises(SolpolpyError, match="degenerate"):
        solve_three_polarizer_brightness(brightness, duplicate_angles, MZP_ANGLES)
