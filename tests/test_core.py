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
