# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from ndcube import NDCollection, NDCube
from pytest import fixture

from solpolpy.constants import VALID_KINDS
from solpolpy.core import _determine_image_shape, add_alpha, determine_input_kind, get_transform_path, resolve
from solpolpy.errors import UnsupportedTransformationError

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = 'WAVE', 'HPLT-TAN', 'HPLN-TAN'
wcs.cdelt = 0.2, 0.5, 0.4
wcs.cunit = 'Angstrom', 'deg', 'deg'
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
@fixture
def npol_ones():
    data_out = []
    data_out.append(("B60.0", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 60, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("B0.0", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 0, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("B-60.0", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': -60, 'OBSRVTRY': 'LASCO'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def fourpol_ones():
    data_out = []
    data_out.append(("B0", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("B45", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 45})))
    data_out.append(("B90", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 90})))
    data_out.append(("B135", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 135})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def mzp_ones():
    data_out = []
    data_out.append(("Bp", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bm", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': -60})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def mzp_data():
    data_out = []
    data_out.append(("Bp", NDCube(np.random.random([50,50]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.random.random([50,50]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bm", NDCube(np.random.random([50,50]), wcs=wcs, meta={'POLAR': -60})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bpb_data():
    data_out = []
    data_out.append(("B", NDCube(np.random.random([50,50]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.random.random([50,50]), wcs=wcs, meta={'POLAR': 'pB'})))
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
    data_out.append(("B", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("alpha", NDCube(np.array([[0]]), wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bpb_ones_no_alpha():
    data_out = []
    data_out.append(("B", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'pB'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Br'})))
    data_out.append(("alpha", NDCube(np.array([[0]])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def btbr_ones_no_alpha():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Br'})))
    # data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def stokes_ones():
    data_out = []
    data_out.append(("Bi", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Bi'})))
    data_out.append(("Bq", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Bq'})))
    data_out.append(("Bu", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'Bu'})))
    return NDCollection(data_out, meta={}, aligned_axes="all")

@fixture
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([[1]]), wcs=wcs, meta={'POLAR': 'pBp'})))
    data_out.append(("alpha", NDCube(np.array([[0]])*u.degree, wcs=wcs)))
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


def test_determine_image_shape():
    data_out = []
    data_out.append(("B", NDCube(np.zeros((20, 20)), wcs=wcs, meta={'POLAR': 'B'})))
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert _determine_image_shape(collection) == (20, 20)


def test_add_alpha(bpb_ones_no_alpha):
    data_out = []
    data_out.append(("B", NDCube(np.zeros((10, 10)), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.zeros((10, 10)), wcs=wcs, meta={'POLAR': 'pB'})))
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert "alpha" not in collection
    out = add_alpha(collection)
    assert "alpha" in out


@pytest.mark.parametrize(
    "fixture_name, kind",
    [("npol_ones", "npol"), ("fourpol_ones", "fourpol"), ("mzp_ones", "mzp"), ("mzp_ones_alpha", "mzp"),
     ("bpb_ones", "bpb"), ("bpb_ones_no_alpha", "bpb"), ("btbr_ones", "btbr"), ("btbr_ones_no_alpha", "btbr"),
     ("stokes_ones", "stokes"), ("bp3_ones", "bp3"), ("bp3_ones_no_alpha", "bp3"), ("bthp_ones", "bthp"),]
)
def test_determine_input_kind(fixture_name, kind, request):
    assert determine_input_kind(request.getfixturevalue(fixture_name)) == kind


def test_determine_input_kind_fail(example_fail):
    with pytest.raises(ValueError):
        determine_input_kind(example_fail)


def test_check_all_paths_resolve(request):
    for kind_source in VALID_KINDS.keys():
        for kind_target in VALID_KINDS.keys():
            if kind_source != kind_target:  # don't need to convert if the input and output systems match
                try:
                    get_transform_path(kind_source, kind_target)
                except UnsupportedTransformationError:
                    pass  # skip any time the transform shouldn't be possible
                else:
                    result = resolve(request.getfixturevalue(f"{kind_source}_ones"), kind_target, out_angles=[0])
                    assert isinstance(result, NDCollection)


def test_btbr_to_npol_missing_out_angles(btbr_ones):
    with pytest.raises(ValueError):
        resolve(btbr_ones, "npol")


def test_imax_effect(mzp_data):
    result = resolve(mzp_data, "MZP", imax_effect=True)
    assert isinstance(result, NDCollection)
    for key in result.keys():
        assert np.sum(result[key].data * mzp_data[key].data) != 0


def test_imax_effect_unsupported_transformation_output(mzp_data):
    with pytest.raises(UnsupportedTransformationError):
        result = resolve(mzp_data, "BpB", imax_effect=True)
        assert isinstance(result, NDCollection)


def test_imax_effect_unsupported_transformation_input(bpb_data):
    with pytest.raises(UnsupportedTransformationError):
        result = resolve(bpb_data, "MZP", imax_effect=True)
        assert isinstance(result, NDCollection)
