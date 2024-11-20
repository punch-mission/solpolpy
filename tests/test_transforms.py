import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from ndcube import NDCollection, NDCube
from pytest import fixture

import solpolpy.transforms as transforms
from solpolpy.errors import MissingAlphaError, SolpolpyError
from tests.fixtures import *

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = "WAVE", "HPLT-TAN", "HPLN-TAN"
wcs.cunit = "Angstrom", "deg", "deg"
wcs.cdelt = 0.2, 0.5, 0.4
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
wcs.cname = "wavelength", "HPC lat", "HPC lon"


def test_bpb_mzp_zeros(bpb_zeros):
    actual = transforms.bpb_to_mzpsolar(bpb_zeros)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_mzp_ones(bpb_ones):
    actual = transforms.bpb_to_mzpsolar(bpb_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([3 / 4]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([3 / 4]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_btbr_ones(btbr_bpb_ones):
    actual = transforms.bpb_to_btbr(btbr_bpb_ones)
    expected_data = []
    expected_data.append(("Bt", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Br", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bt"})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Br"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_bpb_ones(btbr_ones):
    actual = transforms.btbr_to_bpb(btbr_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_btbr_mzp_ways(btbr_ones):
    actual_mzp_direct = transforms.btbr_to_mzpsolar(btbr_ones)
    actual_bpb = transforms.btbr_to_bpb(btbr_ones)
    actual_mzp_indirect = transforms.bpb_to_mzpsolar(actual_bpb)
    for k in list(actual_mzp_direct):
        assert np.allclose(actual_mzp_direct[str(k)].data, actual_mzp_indirect[str(k)].data)


def test_mzp_stokes_ones(mzpsolar_ones):
    actual = transforms.mzpsolar_to_stokes(mzpsolar_ones)
    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def stokes_ones():
    data_out = []
    data_out.append(("I", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bi"})))
    data_out.append(("Q", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bq"})))
    data_out.append(("U", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bu"})))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_stokes_mzp_ones(stokes_ones):
    actual = transforms.stokes_to_mzpsolar(stokes_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([(-np.sqrt(3) + 1) / 4]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([(np.sqrt(3) + 1) / 4]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_mzp_bp3_missing_alpha_errors(mzpsolar_ones):
    with pytest.raises(MissingAlphaError):
        transforms.mzpsolar_to_bp3(mzpsolar_ones)


@fixture()
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pBp"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bp3_mzp_ones(bp3_ones):
    actual = transforms.bp3_to_mzpsolar(bp3_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([-0.5]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_ones_mzp():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bt"})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Br"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_mzp_ones(btbr_ones_mzp):
    actual = transforms.btbr_to_mzpsolar(btbr_ones_mzp)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_bp3_bthp_ones(bp3_ones):
    actual = transforms.bp3_to_bthp(bp3_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("theta", NDCube(np.array([5 * np.pi / 8]), wcs=wcs)))
    expected_data.append(("p", NDCube(np.array([np.sqrt(2)]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_btbr_npol_ones(btbr_ones):
    actual = transforms.btbr_to_npol(btbr_ones, [0, 120, 240] * u.degree)
    expected_data = [(str(0 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     (str(120 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     (str(240 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     ("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_mzp_to_npol_custom():
    """M, Z, P = 0, 1, 0 conversion"""
    input_data = NDCollection(
        [("P", NDCube(np.array([[0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
         ("Z", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
         ("M", NDCube(np.array([[0]]), wcs=wcs, meta={"POLAR": -60 * u.degree}))],
        meta={}, aligned_axes="all")
    actual = transforms.mzpsolar_to_npol(input_data, out_angles=[0, 45, 90] * u.degree)
    expected_data = [(str(0 * u.degree), NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                     (str(45 * u.degree), NDCube(np.array([1 / 3]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
                     (str(90 * u.degree), NDCube(np.array([-1 / 3]), wcs=wcs, meta={"POLAR": 90 * u.degree}))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def fourpol_ones():
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [(str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0})),
                (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45})),
                (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90})),
                (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}))]
    return NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")


def test_fourpol_to_stokes_ones(fourpol_ones):
    actual = transforms.fourpol_to_stokes(fourpol_ones)
    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def mzp_ones_instru():
    input_data = NDCollection(
        [("P", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
         ("Z", NDCube(np.array([[0]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
         ("M", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'}))],
        meta={}, aligned_axes="all")
    return input_data


def test_mzp_mzp_ones_instru(mzp_ones_instru):
    actual = transforms.mzpinstru_to_mzpsolar(mzp_ones_instru)
    expected_data = [
        ("M", NDCube(np.array([[1.01995]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
        ("Z", NDCube(np.array([[0.00041]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
        ("P", NDCube(np.array([[0.97965]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'}))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data, atol=1.e-5)
        assert (actual[str(k)].meta["POLARREF"] == expected[str(k)].meta["POLARREF"])

@fixture()
def mzp_ones_solar():
    input_data = NDCollection(
        [("P", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
         ("Z", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
         ("M", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'}))],
        meta={}, aligned_axes="all")
    return input_data

def test_mzp_mzp_ones_solar(mzp_ones_solar):
    actual = transforms.mzpsolar_to_mzpinstru(mzp_ones_solar)
    expected_data = [
        ("M", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
        ("Z", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
        ("P", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'}))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert (actual[str(k)].meta["POLARREF"] == expected[str(k)].meta["POLARREF"])

@fixture()
def mzpsolar_degenerate():
    data_out = [("M", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": -60*u.degree, "POLAROFF": 60, "POLARREF": 'Instrument'})),
                ("Z", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 0*u.degree,  "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("P", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 60*u.degree, "POLAROFF": -60,  "POLARREF": 'Instrument'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")

def test_mzpsolar_degenerate(mzpsolar_degenerate):
    with pytest.raises(SolpolpyError, match="Conversion matrix is degenerate"):
        transforms.mzpinstru_to_mzpsolar(mzpsolar_degenerate)


@fixture()
def npol_degenerate():
    data_out = [("60.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("0.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("-60.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": -60*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")

def test_npol_degenerate(npol_degenerate):
    with pytest.raises(SolpolpyError, match="Conversion matrix is degenerate"):
        transforms.npol_to_mzpsolar(npol_degenerate)


def test_mask_propagation_works_when_none_provided(fourpol_ones):
    actual = transforms.fourpol_to_stokes(fourpol_ones)
    expected = NDCollection(
        [("I", NDCube(np.array([2]), wcs=wcs, mask=None)),
         ("Q", NDCube(np.array([0]), wcs=wcs, mask=None)),
         ("U", NDCube(np.array([0]), wcs=wcs, mask=None))],
        meta={},
        aligned_axes="all")
    for key in list(expected):
        assert np.allclose(actual[key].data, expected[key].data)
        assert actual[key].mask is None


def test_mask_propagation_works_mixed_normal_case():
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0}, mask=np.zeros(1, dtype=bool))),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert np.equal(actual[str(k)].mask, expected[str(k)].mask)


def test_mask_propagation_works_all_false_normal_case():
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0}, mask=np.zeros(1, dtype=bool))),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert np.equal(actual[str(k)].mask, expected[str(k)].mask)


def test_mask_propagation_works_when_not_all_specified(fourpol_ones):
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0})),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = [("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))),
                     ("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))),
                     ("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool)))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert actual[str(k)].mask is None
