import astropy.units as u
import astropy.wcs
import numpy as np
from ndcube import NDCollection, NDCube
from pytest import fixture

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = "WAVE", "HPLT-TAN", "HPLN-TAN"
wcs.cdelt = 0.2, 0.5, 0.4
wcs.cunit = "Angstrom", "deg", "deg"
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1


@fixture()
def npol_ones():
    data_out = [(str(60 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60*u.degree, "OBSRVTRY": "LASCO"})),
                (str(0 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0*u.degree, "OBSRVTRY": "LASCO"})),
                (str(-60 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60*u.degree, "OBSRVTRY": "LASCO"}))]
    # data_out.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def fourpol_ones():
    data_out = [(str(0 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0*u.degree})),
                (str(45 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 45*u.degree})),
                (str(90 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 90*u.degree})),
                (str(135 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 135*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzp_ones():
    data_out = [("P", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("Z", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 0*u.degree})),
                ("M", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": -60*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzp_ones_other_order():
    data_out = [("Z", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0*u.degree})),
                ("P", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("M", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzp_data():
    data_out = [("P", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("Z", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": 0*u.degree})),
                ("M", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": -60*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_data():
    data_out = [("B", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "pB"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzp_ones_alpha():
    data_out = [("P", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("Z", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 0*u.degree})),
                ("M", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": -60*u.degree})),
                ("alpha", NDCube(np.array([0]), wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_zeros():
    data_out = [("B", NDCube(np.array([0]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.array([0]), wcs=wcs, meta={"POLAR": "pB"})),
                ("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_ones():
    data_out = [("B", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "pB"})),
                ("alpha", NDCube(np.array([[0]]), wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_ones_no_alpha():
    data_out = [("B", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "pB"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def btbr_ones():
    data_out = [("Bt", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "Bt"})),
                ("Br", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "Br"})),
                ("alpha", NDCube(np.array([[0]]) * u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def btbr_ones_no_alpha():
    data_out = [("Bt", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "Bt"})),
                ("Br", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "Br"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def stokes_ones():
    data_out = [("I", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "I"})),
                ("Q", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "Q"})),
                ("U", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "U"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bp3_ones():
    data_out = [("B", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "pB"})),
                ("pBp", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": "pBp"})),
                ("alpha", NDCube(np.array([[0]]) * u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bp3_ones_no_alpha():
    data_out = [("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})),
                ("pBp", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pBp"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bthp_ones():
    data_out = [("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})),
                ("theta", NDCube(np.array([1]) * u.degree, wcs=wcs, meta={"POLAR": "Theta"})),
                ("p", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Degree of Polarization"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def example_fail():
    data_out = [("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})),
                ("Bm", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bm"})),
                ("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")
