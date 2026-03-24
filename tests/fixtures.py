import astropy.units as u
import astropy.wcs
from astropy.io import fits
import numpy as np
from ndcube import NDCollection, NDCube
from pytest import fixture

from solpolpy.util import make_empty_distortion_model

# Solar WCS
wcs_sol = astropy.wcs.WCS(naxis=2)
wcs_sol.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
wcs_sol.wcs.cunit = "deg", "deg"
wcs_sol.wcs.cdelt = 0.5, 0.4
wcs_sol.wcs.crpix = 2, 2
wcs_sol.wcs.crval = 0.5, 1
wcs_sol.wcs.cname = "HPC lon", "HPC lat"

# Celestial WCS
wcs_cel = astropy.wcs.WCS(naxis=2)
wcs_cel.wcs.ctype = "RA---TAN", "DEC--TAN"
wcs_cel.wcs.cunit = "deg", "deg"
wcs_cel.wcs.cdelt = -0.5, 0.4
wcs_cel.wcs.crpix = 2.0, 2.0
wcs_cel.wcs.crval = 10.0, 20.0
wcs_cel.wcs.cname = "RA", "DEC"

hdr = fits.Header()
hdr.update(wcs_sol.to_header())
hdr.update(wcs_cel.to_header(key="A"))

wcs = astropy.wcs.WCS(hdr)

wcs_new = astropy.wcs.WCS(naxis=2)
wcs_new.wcs.ctype = "HPLN-TAN", "HPLT-TAN"
wcs_new.wcs.cunit = "deg", "deg"
wcs_new.wcs.cdelt = 0.5, 0.4
wcs_new.wcs.crpix = 2, 2
wcs_new.wcs.crval = 0.5, 1
wcs_new.wcs.cname = "HPC lon", "HPC lat"

@fixture()
def npol_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [
        (str(60 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 60 * u.degree, "OBSRVTRY": "LASCO"})),
        (str(0 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 0 * u.degree, "OBSRVTRY": "LASCO"})),
        (str(-60 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": -60 * u.degree, "OBSRVTRY": "LASCO"}))]
    data_out.append(("alpha", NDCube(np.full((5,5), 0), wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def fourpol_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [(str(0 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 0 * u.degree})),
                (str(45 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 45 * u.degree})),
                (str(90 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 90 * u.degree})),
                (str(135 * u.degree), NDCube(data, wcs=wcs, meta={"POLAR": 135 * u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpsolar_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("P", NDCube(data, wcs=wcs, meta={"POLAR": 60 * u.degree, "POLARREF": 'Solar'})),
                ("Z", NDCube(data, wcs=wcs, meta={"POLAR": 0 * u.degree, "POLARREF": 'Solar'})),
                ("M", NDCube(data, wcs=wcs, meta={"POLAR": -60 * u.degree, "POLARREF": 'Solar'}))]
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
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("Z", NDCube(data, wcs=wcs, meta={"POLAR": 0 * u.degree, "POLARREF": 'Solar'})),
                ("P", NDCube(data, wcs=wcs, meta={"POLAR": 60 * u.degree, "POLARREF": 'Solar'})),
                ("M", NDCube(data, wcs=wcs, meta={"POLAR": -60 * u.degree, "POLARREF": 'Solar'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_data():
    data_out = [("B", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.random.random([50, 50]), wcs=wcs, meta={"POLAR": "pB"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpsolar_ones_alpha():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("P", NDCube(data, wcs=wcs, meta={"POLAR": 60 * u.degree})),
                ("Z", NDCube(data, wcs=wcs, meta={"POLAR": 0 * u.degree})),
                ("M", NDCube(data, wcs=wcs, meta={"POLAR": -60 * u.degree})),
                ("alpha", NDCube(np.full((5,5), 0)* u.degree, wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_zeros():
    data_out = [("B", NDCube(np.full((5,5), 0), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.full((5,5), 0), wcs=wcs, meta={"POLAR": "pB"})),
                ("alpha", NDCube(np.full((5,5), 0) , wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(data, wcs=wcs, meta={"POLAR": "pB"})),
                ("alpha", NDCube(np.full((5,5), 0), wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bpb_ones_no_alpha():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(data, wcs=wcs, meta={"POLAR": "pB"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def btbr_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("Bt", NDCube(data, wcs=wcs, meta={"POLAR": "Bt"})),
                ("Br", NDCube(data, wcs=wcs, meta={"POLAR": "Br"})),
                ("alpha", NDCube(np.full((5,5), 0) , wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def btbr_ones_no_alpha():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("Bt", NDCube(data, wcs=wcs, meta={"POLAR": "Bt"})),
                ("Br", NDCube(data, wcs=wcs, meta={"POLAR": "Br"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def stokes_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("I", NDCube(data, wcs=wcs, meta={"POLAR": "I"})),
                ("Q", NDCube(data, wcs=wcs, meta={"POLAR": "Q"})),
                ("U", NDCube(data, wcs=wcs, meta={"POLAR": "U"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bp3_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(data, wcs=wcs, meta={"POLAR": "pB"})),
                ("pBp", NDCube(data, wcs=wcs, meta={"POLAR": "pBp"})),
                ("alpha", NDCube(np.full((5,5), 0) , wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bp3_ones_no_alpha():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(data, wcs=wcs, meta={"POLAR": "pB"})),
                ("pBp", NDCube(data, wcs=wcs, meta={"POLAR": "pBp"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def bthp_ones():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("theta", NDCube(data * u.degree, wcs=wcs, meta={"POLAR": "Theta"})),
                ("p", NDCube(data, wcs=wcs, meta={"POLAR": "Degree of Polarization"}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def mzpinstru_ones():
    data, _ = np.mgrid[0:5, 0:5]
    input_data = NDCollection(
        [(
            "P",
            NDCube(data, wcs=wcs_new, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
            ("Z",
             NDCube(data, wcs=wcs_new, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
            ("M",
             NDCube(data, wcs=wcs_new,
                    meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'}))],
        meta={}, aligned_axes="all")
    return input_data


dist = make_empty_distortion_model(50, np.random.random([50, 50]))


@fixture()
def mzpinstru_distortion():
    wcs.cpdis1 = dist[0]
    wcs.cpdis2 = dist[1]
    data_out = [("M", NDCube(np.random.random([50, 50]), wcs=wcs_new,
                             meta={"POLAR": -60 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("Z", NDCube(np.random.random([50, 50]), wcs=wcs_new,
                             meta={"POLAR": 0 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'})),
                ("P", NDCube(np.random.random([50, 50]), wcs=wcs_new,
                             meta={"POLAR": 60 * u.degree, "POLAROFF": 0, "POLARREF": 'Instrument'}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")


@fixture()
def example_fail():
    data, _ = np.mgrid[0:5, 0:5]
    data_out = [("B", NDCube(data, wcs=wcs, meta={"POLAR": "B"})),
                ("Bm", NDCube(data, wcs=wcs, meta={"POLAR": "Bm"})),
                ("alpha", NDCube(np.full((5,5), 0) , wcs=wcs))]
    return NDCollection(data_out, meta={}, aligned_axes="all")
