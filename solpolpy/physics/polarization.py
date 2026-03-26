"""Equation-based helpers for linear polarization transforms."""

from __future__ import annotations

import astropy.units as u
import numpy as np

from solpolpy.errors import SolpolpyError

MZP_ANGLES = np.array([-60.0, 0.0, 60.0]) * u.degree


def as_angle(angles, default_unit=u.radian):
    """Return an input as an angular quantity."""
    quantity = u.Quantity(angles)
    if quantity.unit == u.dimensionless_unscaled:
        quantity = quantity * default_unit
    return quantity


def wrap_pm_pi(angle):
    """Wrap an angle into [-pi, pi) or the unit-equivalent interval."""
    quantity = as_angle(angle, default_unit=u.radian).to(u.radian)
    wrapped = (quantity + np.pi * u.radian) % (2 * np.pi * u.radian) - np.pi * u.radian
    return wrapped.to(quantity.unit)


def wrap_linear_polarization(angle):
    """Wrap a linear-polarization axis angle into [-pi/2, pi/2)."""
    quantity = as_angle(angle, default_unit=u.radian).to(u.radian)
    wrapped = (quantity + (np.pi / 2) * u.radian) % (np.pi * u.radian) - (np.pi / 2) * u.radian
    return wrapped.to(quantity.unit)


def angle_difference_radians(target_angle, source_angle):
    """Return a wrapped target-source angle difference in radians."""
    return wrap_pm_pi(as_angle(target_angle, u.radian) - as_angle(source_angle, u.radian)).to_value(u.radian)


def _as_radians(angles):
    return as_angle(angles, u.radian).to_value(u.radian)


def bp3_to_polarizer_brightness(B, pB, pBp, alpha, polarizer_angles):
    """Evaluate polarizer brightness at the requested polarizer angles."""
    alpha_rad = _as_radians(alpha)
    angle_stack = []
    for angle in u.Quantity(polarizer_angles):
        delta = angle_difference_radians(angle, alpha)
        angle_stack.append(0.5 * (B - pB * np.cos(2 * delta) - pBp * np.sin(2 * delta)))
    return np.stack(angle_stack, axis=0)


def bp3_from_polarizer_brightness(brightness_stack, polarizer_angles, alpha):
    """Recover B, pB, and pBp from polarizer brightness measurements."""
    alpha_rad = _as_radians(alpha)
    brightness_stack = np.asarray(brightness_stack)
    B = (2.0 / 3.0) * np.sum(brightness_stack, axis=0)

    cos_terms = []
    sin_terms = []
    for brightness, angle in zip(brightness_stack, u.Quantity(polarizer_angles), strict=False):
        delta = angle_difference_radians(angle, alpha)
        cos_terms.append(brightness * np.cos(2 * delta))
        sin_terms.append(brightness * np.sin(2 * delta))

    pB = (-4.0 / 3.0) * np.sum(cos_terms, axis=0)
    pBp = (-4.0 / 3.0) * np.sum(sin_terms, axis=0)
    return B, pB, pBp


def _three_polarizer_matrix(source_angles, target_angles, reference_angle=0 * u.degree):
    source = u.Quantity(source_angles)
    target = u.Quantity(target_angles)
    reference = _as_radians(reference_angle)

    rows = []
    for source_angle in source:
        row = []
        for target_angle in target:
            delta = angle_difference_radians(source_angle, target_angle) - reference
            row.append((4.0 * np.cos(delta) ** 2 - 1.0) / 3.0)
        rows.append(np.stack(row, axis=0))
    return np.stack(rows, axis=0)


def solve_three_polarizer_brightness(observed_brightness, observed_angles, solved_angles, reference_angle=0 * u.degree):
    """Solve Equation 44-style systems for the requested output basis."""
    matrix = _three_polarizer_matrix(observed_angles, solved_angles, reference_angle=reference_angle)
    rhs = np.moveaxis(np.asarray(observed_brightness), 0, -1)[..., None]

    if matrix.ndim > 2:
        matrix = np.moveaxis(matrix, (0, 1), (-2, -1))

    try:
        solution = np.matmul(np.linalg.inv(matrix), rhs)[..., 0]
    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            raise SolpolpyError("Conversion matrix is degenerate") from err
        raise

    return np.moveaxis(solution, -1, 0)


def project_three_polarizer_brightness(source_brightness, source_angles, target_angles, reference_angle=0 * u.degree):
    """Project three-polarizer brightness measurements onto new polarizer angles."""
    projected = []
    reference = _as_radians(reference_angle)

    for target_angle in u.Quantity(target_angles):
        value = 0.0
        for brightness, source_angle in zip(source_brightness, u.Quantity(source_angles), strict=False):
            delta = angle_difference_radians(target_angle, source_angle) - reference
            value = value + brightness * ((4.0 * np.cos(delta) ** 2 - 1.0) / 3.0)
        projected.append(value)

    return np.stack(projected, axis=0)
