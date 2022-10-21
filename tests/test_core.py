# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest
from pytest import fixture

import astropy.units as u

from solpolpy.core import determine_input_kind, resolve, sanitize_data_dict
from solpolpy.polarizers import mzp_to_bpb, bpb_to_btbr


@fixture
def example_mzp():
    d = {-60*u.degree: np.array([2]), 0*u.degree: np.array([1]), 60*u.degree: np.array([1.5]),
         "alpha": np.array([5])*u.degree}
    return d


@fixture
def mixed_mzp():
    """MZP data that uses both degrees and radians"""
    d = {-60*u.degree: np.array([2]), 0*u.radian: np.array([1]), np.pi/3*u.radian: np.array([1.5]),
         "alpha": np.array([5])*u.degree}
    return d


@fixture
def unitless_mzp():
    d = {-60: np.array([2]), 0: np.array([1]), 60: np.array([1.5]),
         "alpha": np.array([5])}
    return d


@fixture
def incorrect_units_mzp():
    d = {-60 * u.meter: np.array([2]), 0: np.array([1]), 60: np.array([1.5]),
         "alpha": np.array([5])*u.degree}
    return d


def test_sanitize_mzp_deg2deg(example_mzp):
    """ converts mzp from deg 2 deg"""
    sanitized, use_radians = sanitize_data_dict(example_mzp, u.degree)
    assert not use_radians
    assert example_mzp == sanitized


def test_sanitize_mzp_deg2rad(example_mzp):
    """ converts mzp deg 2 rad"""
    sanitized, use_radians = sanitize_data_dict(example_mzp, u.radian)
    assert not use_radians
    assert set(sanitized.keys()) == {np.pi/3*u.radian, 0*u.radian, -np.pi/3*u.radian, 'alpha'}
    assert sanitized['alpha'][0] == np.pi/36 * u.radian


def test_sanitize_mzp_mixed2rad(mixed_mzp):
    """ converts mzp with mixed units to radians"""
    sanitized, use_radians = sanitize_data_dict(mixed_mzp, u.radian)
    assert use_radians
    assert set(sanitized.keys()) == {np.pi / 3 * u.radian, 0 * u.radian, -np.pi / 3 * u.radian, 'alpha'}
    assert sanitized['alpha'][0] == np.pi / 36 * u.radian


def test_sanitize_mzp_mixed2deg(mixed_mzp):
    """ converts mzp mixed units to degrees"""
    sanitized, use_radians = sanitize_data_dict(mixed_mzp, u.degree)
    assert use_radians
    assert -60 * u.degree in sanitized
    assert (0*u.rad) * ((180*u.degree) / (np.pi* u.radian)) in sanitized
    assert (np.pi/3*u.rad) * ((180*u.degree) / (np.pi* u.radian)) in sanitized
    assert sanitized['alpha'][0] == 5 * u.degree


def test_sanitize_mzp_unitless2deg(unitless_mzp):
    """ converts unitless mzp to degrees"""
    sanitized, use_radians = sanitize_data_dict(unitless_mzp, u.degree)
    assert not use_radians
    assert -60 * u.degree in sanitized
    assert 0 * u.degree in sanitized
    assert 60 * u.degree in sanitized
    assert sanitized['alpha'][0] == 5 * u.degree


def test_sanitize_mzp_unitless2rad(unitless_mzp):
    """ converts unitless mzp to radians"""
    sanitized, use_radians = sanitize_data_dict(unitless_mzp, u.radian)
    assert not use_radians
    assert -np.pi/3 * u.radian in sanitized
    assert 0 * u.radian in sanitized
    assert np.pi/3 * u.radian in sanitized
    assert sanitized['alpha'][0] == np.pi/36 * u.radian


def test_sanitize_wrong_unit_request(example_mzp):
    """ tries to convert to a nonsensical unit"""
    with pytest.raises(RuntimeError):
        sanitize_data_dict(example_mzp, 'fizbop')


def test_sanitize_wrong_unit_mzp(incorrect_units_mzp):
    """ key has the wrong unit (meters instead of an angle unit)"""
    with pytest.raises(RuntimeError):
        sanitize_data_dict(incorrect_units_mzp, u.degree)


def test_sanitize_wrong_unit_alpha():
    """ alpha has non-angle unit of meter"""
    d = {"alpha": np.array([5])*u.meter}
    with pytest.raises(RuntimeError):
        sanitize_data_dict(d, u.degree)


def test_sanitize_nonnumeric_alpha():
    """ alpha is nonnumeric so can't sanitize"""
    d = {"alpha": np.array(["fail"])}
    with pytest.raises(RuntimeError):
        sanitize_data_dict(d, u.degree)


def test_sanitize_radian2degree_alpha():
    """ alpha gets converted from radians to degrees"""
    d = {"alpha": np.array([0])*u.radian}
    sanitized, use_radians = sanitize_data_dict(d, u.degree)
    assert use_radians
    assert sanitized['alpha'][0] == 0 * u.degree


def test_sanitize_radian2radian_alpha():
    """ alpha gets converted from radians to radians"""
    d = {"alpha": np.array([0])*u.radian}
    sanitized, use_radians = sanitize_data_dict(d, u.radian)
    assert use_radians
    assert sanitized['alpha'][0] == 0 * u.radian


def test_determine_input_kind_mzp(example_mzp):
    assert determine_input_kind(example_mzp), "MZP"


def test_determine_input_kind():
    d = {"B": np.array([0]), "pB": np.array([1])}
    determined_kind = determine_input_kind(d)
    assert determined_kind, "BpB"


def test_determine_input_kind_fail():
    d = {"B": np.array([0]), "M": np.array([1])}
    with pytest.raises(ValueError):
        determine_input_kind(d)


def test_bpb_to_btbr():
    d = {"B": np.array([2]), "pB": np.array([1])}
    result = resolve(d, "BtBr")
    assert isinstance(result, dict)


def test_mzp_to_btbr(example_mzp):
    result = resolve(example_mzp, "BtBr")
    assert isinstance(result, dict)


def test_STEREO_FITS():
    """ingest STEREO data nd test does something"""
    file_list=["./stereo1.fts","./stereo2.fts","./stereo3.fts"]
    result = resolve(file_list, "BtBr")
    assert isinstance(result, dict)