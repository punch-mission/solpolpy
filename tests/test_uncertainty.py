import astropy.units as u
import astropy.wcs
import numpy as np
from ndcube import NDCollection, NDCube
from pytest import fixture

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = 'WAVE', 'HPLT-TAN', 'HPLN-TAN'
wcs.cunit = 'Angstrom', 'deg', 'deg'
wcs.cdelt = 0.2, 0.5, 0.4
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
wcs.cname = 'wavelength', 'HPC lat', 'HPC lon'

@fixture
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'B'})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pB'})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={'POLAR': 'pBp'})))
    data_out.append(("alpha", NDCube(np.array([0])*u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")
