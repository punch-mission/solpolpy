import astropy.units as u
import astropy.wcs
import numpy as np
from ndcube import NDCollection, NDCube
from pytest import fixture

from solpolpy.util import make_empty_distortion_model

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = "WAVE", "HPLT-TAN", "HPLN-TAN"
wcs.cdelt = 0.2, 0.5, 0.4
wcs.cunit = "Angstrom", "deg", "deg"
wcs.crpix = 2, 2, 2
wcs.crval = 0, 0, 0


@fixture()
def npol_ones():
    data_out = [
        (str(60 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "OBSRVTRY": "LASCO"})),
        (str(0 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "OBSRVTRY": "LASCO"})),
        (str(-60 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "OBSRVTRY": "LASCO"}))]
    # data_out.append(("alpha", NDCube(np.array([0])*u.radian, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def fourpol_ones():
    data_out = [(str(0 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                (str(45 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
                (str(90 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 90 * u.degree})),
                (str(135 * u.degree), NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 135 * u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpsolar_ones():
    data_out = [("P", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLARREF": 'Solar'})),
                ("Z", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLARREF": 'Solar'})),
                ("M", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLARREF": 'Solar'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpsolar_degenerate():
    data_out = [("P", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs,
                             meta={"POLAR": 60 * u.degree, "POLAROFF": -60, "POLARREF": 'Instrument'})),
                ("Z", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs,
                             meta={"POLAR": 0 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("M", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs,
                             meta={"POLAR": -60 * u.degree, "POLAROFF": 60, "POLARREF": 'Instrument'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzp_ones_other_order():
    data_out = [("Z", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLARREF": 'Solar'})),
                ("P", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLARREF": 'Solar'})),
                ("M", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLARREF": 'Solar'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_data():
    data_out = [("B", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "pB"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpsolar_ones_alpha():
    data_out = [("P", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
                ("Z", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                ("M", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
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
def mzpinstru_ones():
    input_data = NDCollection(
        [(
            "P",
            NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
            ("Z",
             NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
            ("M",
             NDCube(np.array([[1]]), wcs=wcs,
                    meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'}))],
        meta={}, aligned_axes="all")
    return input_data


dist = make_empty_distortion_model(50, np.random.random([50, 50]))


@fixture()
def mzpinstru_distortion():
    wcs.cpdis1 = dist[0]
    wcs.cpdis2 = dist[1]
    data_out = [("M", NDCube(np.random.random([50, 50]), wcs=wcs,
                             meta={"POLAR": -60 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("Z", NDCube(np.random.random([50, 50]), wcs=wcs,
                             meta={"POLAR": 0 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("P", NDCube(np.random.random([50, 50]), wcs=wcs,
                             meta={"POLAR": 60 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def example_fail():
    data_out = [("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})),
                ("Bm", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bm"})),
                ("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")
