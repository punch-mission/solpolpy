# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest
from pytest import fixture
import astropy.wcs
from ndcube import NDCollection, NDCube
import os

import astropy.units as u

from solpolpy.core import (determine_input_kind,
                          resolve,
                          )
# from solpolpy.polarizers import mzp_to_bpb, bpb_to_btbr


wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = 'WAVE', 'HPLT-TAN', 'HPLN-TAN'
wcs.cdelt = 0.2, 0.5, 0.4
wcs.cunit = 'Angstrom', 'deg', 'deg'
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
@fixture
def npol_ones():
    data_out = []
    data_out.append(("angle_1", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("angle_2", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("angle_3", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60, 'OBSRVTRY': 'LASCO'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def four_pol_ones():
    data_out = []
    data_out.append(("angle_1", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("angle_2", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 45})))
    data_out.append(("angle_3", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 90})))
    data_out.append(("angle_4", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 135})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def mzp_ones():
    data_out = []
    data_out.append(("Bp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60})))
    return NDCollection(data_out, meta={}, aligned_axes="all")
@fixture
def mzp_ones_alpha():
    data_out = []
    data_out.append(("Bp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60})))
    data_out.append(("alpha", NDCube(np.array([0]), wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bpb_ones_no_alpha():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Br'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def btbr_ones_no_alpha():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Br'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def stokes_ones():
    data_out = []
    data_out.append(("Bi", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bi'})))
    data_out.append(("Bq", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bq'})))
    data_out.append(("Bu", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bu'})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pBp'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bp3_ones_no_alpha():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pBp'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bthp_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("theta", NDCube(np.array([1])*u.degree, wcs=wcs, meta={'POLAR': 'Theta'})))
    data_out.append(("p", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Degree of Polarization'})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def example_fail():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bm'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


'''
    Tests to determine_input_kind
'''


def test_determine_input_kind_npol(npol_ones):
    assert determine_input_kind(npol_ones), "npol"

def test_determine_input_kind_fourpol(four_pol_ones):
    assert determine_input_kind(four_pol_ones), "fourpol"

def test_determine_input_kind_mzp(mzp_ones):
    assert determine_input_kind(mzp_ones), "MZP"

def test_determine_input_kind_mzp_alpha(mzp_ones_alpha):
    assert determine_input_kind(mzp_ones_alpha), "MZP"

def test_determine_input_kind_bpb(bpb_ones):
    assert determine_input_kind(bpb_ones), "BpB"

def test_determine_input_kind_bpb_noa(bpb_ones_no_alpha):
    assert determine_input_kind(bpb_ones_no_alpha), "BpB"

def test_determine_input_kind_btbr(btbr_ones):
    assert determine_input_kind(btbr_ones), "BtBr"

def test_determine_input_kind_btbr_noa(btbr_ones_no_alpha):
    assert determine_input_kind(btbr_ones_no_alpha), "BtBr"

def test_determine_input_kind_stokes(stokes_ones):
    assert determine_input_kind(stokes_ones), "Stokes"

def test_determine_input_kind_bp3(bp3_ones):
    assert determine_input_kind(bp3_ones), "Bp3"

def test_determine_input_kind_bp3_noa(bp3_ones_no_alpha):
    assert determine_input_kind(bp3_ones_no_alpha), "Bp3"

def test_determine_input_kind_bthp(bthp_ones):
    assert determine_input_kind(bthp_ones), "Bthp"

def test_determine_input_kind_fail(example_fail):
    with pytest.raises(ValueError):
        determine_input_kind(example_fail)


'''
    Tests to check the resolve
'''

def test_npol_to_mzp(npol_ones):
    result = resolve(npol_ones, "MZP")
    assert isinstance(result, NDCollection)

def test_mzp_to_bpb(mzp_ones_alpha):
    result = resolve(mzp_ones_alpha, "BpB")
    assert isinstance(result, NDCollection)

def test_bpb_to_mzp(bpb_ones):
    result = resolve(bpb_ones, "MZP")
    assert isinstance(result, NDCollection)

def test_bpb_to_btbr(bpb_ones):
    result = resolve(bpb_ones, "BtBr")
    assert isinstance(result, NDCollection)

def test_btbr_to_bpb(btbr_ones):
    result = resolve(btbr_ones, "BpB")
    assert isinstance(result, NDCollection)

def test_mzp_to_stokes(mzp_ones):
    result = resolve(mzp_ones, "Stokes")
    assert isinstance(result, NDCollection)

def test_stokes_to_mzp(stokes_ones):
    result = resolve(stokes_ones, "MZP")
    assert isinstance(result, NDCollection)

def test_mzp_to_bp3(mzp_ones_alpha):
    result = resolve(mzp_ones_alpha, "Bp3")
    assert isinstance(result, NDCollection)

def test_bp3_to_mzp(bp3_ones):
    result = resolve(bp3_ones, "MZP")
    assert isinstance(result, NDCollection)

def test_btbr_to_mzp(btbr_ones):
    result = resolve(btbr_ones, "MZP")
    assert isinstance(result, NDCollection)

def test_bp3_to_bthp(bp3_ones):
    result = resolve(bp3_ones, "Bthp")
    assert isinstance(result, NDCollection)


# def test_determine_input_kind():
#     d = {"B": np.array([0]), "pB": np.array([1])}
#     determined_kind = determine_input_kind(d)
#     assert determined_kind, "BpB"
#
#
# def test_determine_input_kind_fail():
#     d = {"B": np.array([0]), "M": np.array([1])}
#     with pytest.raises(ValueError):
#         determine_input_kind(d)
# @fixture
# def example_mzp():
#     d = {-60*u.degree: np.array([2]), 0*u.degree: np.array([1]), 60*u.degree: np.array([1.5]),
#          "alpha": np.array([5])*u.degree}
#     return d
#
#
# @fixture
# def mixed_mzp():
#     """MZP data that uses both degrees and radians"""
#     d = {-60*u.degree: np.array([2]), 0*u.radian: np.array([1]), np.pi/3*u.radian: np.array([1.5]),
#          "alpha": np.array([5])*u.degree}
#     return d
#
#
# @fixture
# def unitless_mzp():
#     d = {-60: np.array([2]), 0: np.array([1]), 60: np.array([1.5]),
#          "alpha": np.array([5])}
#     return d
#
#
# @fixture
# def incorrect_units_mzp():
#     d = {-60 * u.meter: np.array([2]), 0: np.array([1]), 60: np.array([1.5]),
#          "alpha": np.array([5])*u.degree}
#     return d
#
#
# def test_sanitize_mzp_deg2deg(example_mzp):
#     """ converts mzp from deg 2 deg"""
#     sanitized, use_radians = sanitize_data_dict(example_mzp, u.degree)
#     assert not use_radians
#     assert example_mzp == sanitized
#
#
# def test_sanitize_mzp_deg2rad(example_mzp):
#     """ converts mzp deg 2 rad"""
#     sanitized, use_radians = sanitize_data_dict(example_mzp, u.radian)
#     assert not use_radians
#     assert set(sanitized.keys()) == {np.pi/3*u.radian, 0*u.radian, -np.pi/3*u.radian, 'alpha'}
#     assert sanitized['alpha'][0] == np.pi/36 * u.radian
#
#
# def test_sanitize_mzp_mixed2rad(mixed_mzp):
#     """ converts mzp with mixed units to radians"""
#     sanitized, use_radians = sanitize_data_dict(mixed_mzp, u.radian)
#     assert use_radians
#     assert set(sanitized.keys()) == {np.pi / 3 * u.radian, 0 * u.radian, -np.pi / 3 * u.radian, 'alpha'}
#     assert sanitized['alpha'][0] == np.pi / 36 * u.radian
#
#
# def test_sanitize_mzp_mixed2deg(mixed_mzp):
#     """ converts mzp mixed units to degrees"""
#     sanitized, use_radians = sanitize_data_dict(mixed_mzp, u.degree)
#     assert use_radians
#     assert -60 * u.degree in sanitized
#     assert (0*u.rad) * ((180*u.degree) / (np.pi* u.radian)) in sanitized
#     assert (np.pi/3*u.rad) * ((180*u.degree) / (np.pi* u.radian)) in sanitized
#     assert sanitized['alpha'][0] == 5 * u.degree
#
#
# def test_sanitize_mzp_unitless2deg(unitless_mzp):
#     """ converts unitless mzp to degrees"""
#     sanitized, use_radians = sanitize_data_dict(unitless_mzp, u.degree)
#     assert not use_radians
#     assert -60 * u.degree in sanitized
#     assert 0 * u.degree in sanitized
#     assert 60 * u.degree in sanitized
#     assert sanitized['alpha'][0] == 5 * u.degree
#
#
# def test_sanitize_mzp_unitless2rad(unitless_mzp):
#     """ converts unitless mzp to radians"""
#     sanitized, use_radians = sanitize_data_dict(unitless_mzp, u.radian)
#     assert not use_radians
#     assert -np.pi/3 * u.radian in sanitized
#     assert 0 * u.radian in sanitized
#     assert np.pi/3 * u.radian in sanitized
#     assert sanitized['alpha'][0] == np.pi/36 * u.radian
#
#
# def test_sanitize_wrong_unit_request(example_mzp):
#     """ tries to convert to a nonsensical unit"""
#     with pytest.raises(RuntimeError):
#         sanitize_data_dict(example_mzp, 'fizbop')
#
#
# def test_sanitize_wrong_unit_mzp(incorrect_units_mzp):
#     """ key has the wrong unit (meters instead of an angle unit)"""
#     with pytest.raises(RuntimeError):
#         sanitize_data_dict(incorrect_units_mzp, u.degree)
#
#
# def test_sanitize_wrong_unit_alpha():
#     """ alpha has non-angle unit of meter"""
#     d = {"alpha": np.array([5])*u.meter}
#     with pytest.raises(RuntimeError):
#         sanitize_data_dict(d, u.degree)
#
#
# def test_sanitize_nonnumeric_alpha():
#     """ alpha is nonnumeric so can't sanitize"""
#     d = {"alpha": np.array(["fail"])}
#     with pytest.raises(RuntimeError):
#         sanitize_data_dict(d, u.degree)
#
#
# def test_sanitize_radian2degree_alpha():
#     """ alpha gets converted from radians to degrees"""
#     d = {"alpha": np.array([0])*u.radian}
#     sanitized, use_radians = sanitize_data_dict(d, u.degree)
#     assert use_radians
#     assert sanitized['alpha'][0] == 0 * u.degree
#
#
# def test_sanitize_radian2radian_alpha():
#     """ alpha gets converted from radians to radians"""
#     d = {"alpha": np.array([0])*u.radian}
#     sanitized, use_radians = sanitize_data_dict(d, u.radian)
#     assert use_radians
#     assert sanitized['alpha'][0] == 0 * u.radian
#
#
# def test_determine_input_kind_mzp(example_mzp):
#     assert determine_input_kind(example_mzp), "MZP"
#
#
# def test_determine_input_kind():
#     d = {"B": np.array([0]), "pB": np.array([1])}
#     determined_kind = determine_input_kind(d)
#     assert determined_kind, "BpB"
#
#
# def test_determine_input_kind_fail():
#     d = {"B": np.array([0]), "M": np.array([1])}
#     with pytest.raises(ValueError):
#         determine_input_kind(d)
#
#
# def test_bpb_to_btbr():
#     d = {"B": np.array([2]), "pB": np.array([1])}
#     result = resolve(d, "BtBr")
#     assert isinstance(result, dict)
#
#
# def test_mzp_to_btbr(example_mzp):
#     result = resolve(example_mzp, "BtBr")
#     assert isinstance(result, dict)
#
#
# def test_STEREO_triplet():
#     """ingest STEREO data nd test does something"""
#     TESTDATA_DIR = os.path.dirname(__file__)
#     path_to_test_files=TESTDATA_DIR+'/test_support_files/'
#     file_list=[path_to_test_files+"stereo_0.fts",
#                path_to_test_files+"stereo_120.fts",
#                path_to_test_files+"stereo_240.fts"]
#     result = convert_image_list_to_dict(file_list)
#     assert isinstance(result, dict)
#     assert 0*u.degree in result
#     assert 120*u.degree in result
#     assert 240*u.degree in result
#
#
# def test_LASCO_triplet():
#     """ingest STEREO data nd test does something"""
#     TESTDATA_DIR = os.path.dirname(__file__)
#     path_to_test_files=TESTDATA_DIR+'/test_support_files/'
#     file_list=[path_to_test_files+"lasco_-60.fts",
#                path_to_test_files+"lasco_+60.fts",
#                path_to_test_files+"lasco_0.fts"]
#     result = convert_image_list_to_dict(file_list)
#     assert isinstance(result, dict)
#     assert 0*u.degree in result
#     assert -60*u.degree in result
#     assert 60*u.degree in result
#
#
# def test_LASCO_BpB():
#     """ingest STEREO data nd test does something"""
#     TESTDATA_DIR = os.path.dirname(__file__)
#     path_to_test_files=TESTDATA_DIR+'/test_support_files/'
#     file_list=[path_to_test_files+"lasco_-60.fts",
#                path_to_test_files+"lasco_clear.fts"]
#     result = convert_image_list_to_dict(file_list)
#     assert isinstance(result, dict)
#     assert 'pB' in result
#     assert 'B' in result
#     assert 'alpha' in result
#
#
# @pytest.mark.parametrize("B, B_theta, theta, alpha, expected",
#                             [(0, 0.5, 0*u.degree, 0*u.degree, -1),
#                             (0, 0.5, 0*u.radian, 0*u.radian, -1),
#                             (0, 0.5, 0*u.radian, 0*u.degree, -1),
#                             (0, 0.5, 1*u.degree, 1*u.degree, -1),
#                             (0, 0.5, 1*u.degree, 2*u.degree, -1.00061)])
#
# def test_pB_from_single_angle_function(B, B_theta, theta, alpha, expected):
#     """test pB from single angle fn in core"""
#
#     output=pB_from_single_angle(B, B_theta, theta, alpha)
#     np.testing.assert_allclose(output, expected, rtol=1e-05)
