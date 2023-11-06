import numpy as np
import astropy.units as u
import astropy.wcs
from ndcube import NDCollection, NDCube
import pytest
from pytest import fixture

import solpolpy.polarizers as pol

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = 'WAVE', 'HPLT-TAN', 'HPLN-TAN'
wcs.cunit = 'Angstrom', 'deg', 'deg'
wcs.cdelt = 0.2, 0.5, 0.4
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
wcs.cname = 'wavelength', 'HPC lat', 'HPC lon'


@fixture
def npol_mzp_zeros():
    data_out = []
    data_out.append(("angle_1", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 60.0, 'OBSRVTRY': 'STEREO_A'})))
    data_out.append(("angle_2", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 0.0, 'OBSRVTRY': 'STEREO_A'})))
    data_out.append(("angle_3", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': -60.0, 'OBSRVTRY': 'STEREO_A'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_npol_mzp_zeros(npol_mzp_zeros):
    actual = pol.npol_to_mzp(npol_mzp_zeros)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def npol_mzp_ones():
    data_out = []
    data_out.append(("angle_1", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60, 'OBSRVTRY': 'STEREO_B'})))
    data_out.append(("angle_2", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0, 'OBSRVTRY': 'STEREO_B'})))
    data_out.append(("angle_3", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60, 'OBSRVTRY': 'STEREO_B'})))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_npol_mzp_ones(npol_mzp_ones):
    actual = pol.npol_to_mzp(npol_mzp_ones)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def mzp_ones_alpha():
    data_out = []
    data_out.append(("angle_1", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("angle_2", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("angle_3", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60, 'OBSRVTRY': 'LASCO'})))
    data_out.append(("alpha", NDCube(np.array([0]), wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_npol_mzp_ones_alpha(mzp_ones_alpha):
    actual = pol.npol_to_mzp(mzp_ones_alpha)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def mzp_zeros():
    data_out = []
    data_out.append(("Bm", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bp", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': -60})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_mzp_bpb_zeros(mzp_zeros):
    actual = pol.mzp_to_bpb(mzp_zeros)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def mzp_ones():
    data_out = []
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_mzp_bpb_ones(mzp_ones):
    actual = pol.mzp_to_bpb(mzp_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.zeros(1), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def bpb_zeros():
    data_out = []
    data_out.append(("B", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([0]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_mzp_zeros(bpb_zeros):
    actual = pol.bpb_to_mzp(bpb_zeros)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_mzp_ones(bpb_ones):
    actual = pol.bpb_to_mzp(bpb_ones)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([3/4]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([3/4]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def btbr_bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_btbr_ones(btbr_bpb_ones):
    actual = pol.bpb_to_btbr(btbr_bpb_ones)
    expected_data = []
    expected_data.append(("Bt", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Br", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Br'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_bpb_ones(btbr_ones):
    actual = pol.btbr_to_bpb(btbr_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def mzp_ones():
    data_out = []
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_mzp_stokes_ones(mzp_ones):
    actual = pol.mzp_to_stokes(mzp_ones)
    expected_data = []
    expected_data.append(("Bi", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("Bq", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bu", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def stokes_ones():
    data_out = []
    data_out.append(("Bi", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bi'})))
    data_out.append(("Bq", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bq'})))
    data_out.append(("Bu", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bu'})))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_stokes_mzp_ones(stokes_ones):
    actual = pol.stokes_to_mzp(stokes_ones)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([(-np.sqrt(3)+1)/4]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([(np.sqrt(3)+1)/4]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def mzp_ones():
    data_out = []
    data_out.append(("Bm", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 60})))
    data_out.append(("Bz", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("Bp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': -60})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_mzp_bp3_ones(mzp_ones):
    actual = pol.mzp_to_bp3(mzp_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.zeros(1), wcs=wcs)))
    expected_data.append(("pBp", NDCube(np.zeros(1), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pBp'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bp3_mzp_ones(bp3_ones):
    actual = pol.bp3_to_mzp(bp3_ones)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([-0.5]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def btbr_ones_mzp():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Br'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_mzp_ones(btbr_ones_mzp):
    actual = pol.btbr_to_mzp(btbr_ones_mzp)
    expected_data = []
    expected_data.append(("Bm", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bz", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Bp", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pBp'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bp3_bthp_ones(bp3_ones):
    actual = pol.bp3_to_bthp(bp3_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("theta", NDCube(np.array([5*np.pi/8]), wcs=wcs)))
    expected_data.append(("p", NDCube(np.array([np.sqrt(2)]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Bt'})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'Br'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_npol_ones(btbr_ones):
    actual = pol.btbr_to_npol(btbr_ones,[0,120,240])
    expected_data = []
    expected_data.append(("B0", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("B120", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("B240", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture
def fourpol_ones():
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = 'ONE'
    wcs.cunit = 'deg'
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = 'ONE'

    data_out = []
    data_out.append(("B0", NDCube(data=np.array([1]), wcs=wcs, meta={'POLAR': 0})))
    data_out.append(("B45", NDCube(data=np.array([1]), wcs=wcs, meta={'POLAR': 45})))
    data_out.append(("B90", NDCube(data=np.array([1]), wcs=wcs, meta={'POLAR': 90})))
    data_out.append(("B135", NDCube(data=np.array([1]), wcs=wcs, meta={'POLAR': 135})))
    return NDCollection(key_data_pairs=data_out, meta={}, aligned_axes='all')


def test_fourpol_to_stokes_ones(fourpol_ones):
    actual = pol.fourpol_to_stokes(fourpol_ones)
    expected_data = []
    expected_data.append(("Bi", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("Bq", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Bu", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
