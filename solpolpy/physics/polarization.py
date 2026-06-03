"""Equation-based helpers for linear polarization transforms."""

from __future__ import annotations

import astropy.units as u
import numpy as np

from solpolpy.errors import SolpolpyError

MZP_ANGLES = np.array([-60.0, 0.0, 60.0]) * u.degree


def as_angle(angles, default_unit: u.Unit = u.radian) -> u.Quantity:
    """Return an input as an angular quantity.

    Scalars and arrays without units are interpreted in ``default_unit``. This
    helper is intentionally more permissive than ``quantity_input`` because
    FITS metadata and older call sites can provide plain numeric angle values.
    """
    if not isinstance(angles, u.Quantity) and hasattr(angles, "value"):
        angles = angles.value
    quantity = u.Quantity(angles)
    if quantity.unit == u.dimensionless_unscaled:
        quantity = quantity * default_unit
    return quantity


@u.quantity_input(angle=u.radian)
def wrap_pm_pi(angle: u.Quantity):
    """Wrap an angle into [-pi, pi) or the unit-equivalent interval."""
    input_unit = angle.unit
    quantity = angle.to(u.radian)
    wrapped = (quantity + np.pi * u.radian) % (2 * np.pi * u.radian) - np.pi * u.radian
    return wrapped.to(input_unit)


@u.quantity_input(angle=u.radian)
def wrap_linear_polarization(angle: u.Quantity):
    """Wrap a linear-polarization axis angle into [-pi/2, pi/2)."""
    input_unit = angle.unit
    quantity = angle.to(u.radian)
    wrapped = (quantity + (np.pi / 2) * u.radian) % (np.pi * u.radian) - (np.pi / 2) * u.radian
    return wrapped.to(input_unit)


@u.quantity_input(target_angle=u.radian, source_angle=u.radian)
def calculate_angle_difference(target_angle: u.Quantity, source_angle: u.Quantity):
    """Return the wrapped ``target_angle - source_angle`` difference."""
    return wrap_pm_pi(target_angle - source_angle)


def _as_radians(angles: u.Quantity) -> np.ndarray:
    """Return angle values in radians for numeric trigonometric operations."""
    return as_angle(angles, u.radian).to_value(u.radian)


def _angle_difference_radians(target_angle: u.Quantity, source_angle: u.Quantity) -> np.ndarray:
    """Return wrapped angle-difference values in radians."""
    return calculate_angle_difference(as_angle(target_angle, u.radian), as_angle(source_angle, u.radian)).to_value(u.radian)


def bp3_to_polarizer_brightness(
    B: np.ndarray,
    pB: np.ndarray,
    pBp: np.ndarray,
    alpha: u.Quantity,
    polarizer_angles: u.Quantity,
) -> np.ndarray:
    """Evaluate polarizer brightness at the requested polarizer angles."""
    angle_stack = []
    for angle in u.Quantity(polarizer_angles):
        delta = _angle_difference_radians(angle, alpha)
        angle_stack.append(0.5 * (B - pB * np.cos(2 * delta) - pBp * np.sin(2 * delta)))
    return np.stack(angle_stack, axis=0)


def bp3_from_polarizer_brightness(
    brightness_stack: np.ndarray,
    polarizer_angles: u.Quantity,
    alpha: u.Quantity,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Recover B, pB, and pBp from polarizer brightness measurements."""
    brightness_stack = np.asarray(brightness_stack)
    B = (2.0 / 3.0) * np.sum(brightness_stack, axis=0)

    cos_terms = []
    sin_terms = []
    for brightness, angle in zip(brightness_stack, u.Quantity(polarizer_angles), strict=False):
        delta = _angle_difference_radians(angle, alpha)
        cos_terms.append(brightness * np.cos(2 * delta))
        sin_terms.append(brightness * np.sin(2 * delta))

    pB = (-4.0 / 3.0) * np.sum(cos_terms, axis=0)
    pBp = (-4.0 / 3.0) * np.sum(sin_terms, axis=0)
    return B, pB, pBp


def _three_polarizer_matrix(
    source_angles: u.Quantity,
    target_angles: u.Quantity,
    reference_angle: u.Quantity = 0 * u.degree,
) -> np.ndarray:
    """Return the Equation 44 projection matrix for three polarizer angles."""
    source = u.Quantity(source_angles)
    target = u.Quantity(target_angles)
    reference = _as_radians(reference_angle)

    rows = []
    for source_angle in source:
        row = []
        for target_angle in target:
            delta = _angle_difference_radians(source_angle, target_angle) - reference
            row.append((4.0 * np.cos(delta) ** 2 - 1.0) / 3.0)
        rows.append(np.stack(row, axis=0))
    return np.stack(rows, axis=0)


def solve_to_three_polarizer_brightness(
    observed_brightness: np.ndarray,
    observed_angles: u.Quantity,
    solved_angles: u.Quantity,
    reference_angle: u.Quantity = 0 * u.degree,
) -> np.ndarray:
    """
    Solve Equation 44-style systems for the requested output basis using matrix inversion.
    Estimating B_i on RHS of this equation.
    """
    matrix = _three_polarizer_matrix(observed_angles, solved_angles, reference_angle=reference_angle)
    rhs = np.moveaxis(np.asarray(observed_brightness), 0, -1)[..., None]

    if matrix.ndim > 2:
        # Pixelwise reference angles produce a matrix per pixel. Move the
        # matrix axes behind the image axes so numpy can solve all pixels at
        # once with batched matrix multiplication.
        matrix = np.moveaxis(matrix, (0, 1), (-2, -1))

    try:
        solution = np.matmul(np.linalg.inv(matrix), rhs)[..., 0]
    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            raise SolpolpyError("Conversion matrix is degenerate") from err
        raise

    return np.moveaxis(solution, -1, 0)


def solve_from_three_polarizer_brightness(
    source_brightness: np.ndarray,
    source_angles: u.Quantity,
    target_angles: u.Quantity,
    reference_angle: u.Quantity = 0 * u.degree,
) -> np.ndarray:
    """
    Solve Equation 45-style systems brightness measurements onto new polarizer angles.
    Estimating B_phi on LHS of this equation.
    """
    projected = []
    reference = _as_radians(reference_angle)

    for target_angle in u.Quantity(target_angles):
        value = 0.0
        for brightness, source_angle in zip(source_brightness, u.Quantity(source_angles), strict=False):
            delta = _angle_difference_radians(target_angle, source_angle) - reference
            value = value + brightness * ((4.0 * np.cos(delta) ** 2 - 1.0) / 3.0)
        projected.append(value)

    return np.stack(projected, axis=0)
