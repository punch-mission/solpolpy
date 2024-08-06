"""pytest test suite for the polarizers module of solpolpy
"""

import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from ndcube import NDCollection, NDCube

from solpolpy.core import _determine_image_shape, add_alpha, determine_input_kind, get_transform_path, resolve
from solpolpy.errors import UnsupportedTransformationError
from solpolpy.transforms import System
from tests.fixtures import *

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = "WAVE", "HPLT-TAN", "HPLN-TAN"
wcs.cdelt = 0.2, 0.5, 0.4
wcs.cunit = "Angstrom", "deg", "deg"
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1


def test_determine_image_shape():
    data_out = [("B", NDCube(np.zeros((20, 20)), wcs=wcs, meta={"POLAR": "B"}))]
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert _determine_image_shape(collection) == (20, 20)


def test_add_alpha(bpb_ones_no_alpha):
    data_out = [("B", NDCube(np.zeros((10, 10)), wcs=wcs, meta={"POLAR": "B"})),
                ("pB", NDCube(np.zeros((10, 10)), wcs=wcs, meta={"POLAR": "pB"}))]
    collection = NDCollection(data_out, meta={}, aligned_axes="all")
    assert "alpha" not in collection
    out = add_alpha(collection)
    assert "alpha" in out


@pytest.mark.parametrize(
    "fixture_name, kind",
    [("npol_ones", System.npol),
     ("fourpol_ones", System.fourpol),
     ("mzp_ones", System.mzp),
     ("mzp_ones_alpha", System.mzp),
     ("bpb_ones", System.bpb),
     ("bpb_ones_no_alpha", System.bpb),
     ("btbr_ones", System.btbr),
     ("btbr_ones_no_alpha", System.btbr),
     ("stokes_ones", System.stokes),
     ("bp3_ones", System.bp3),
     ("bp3_ones_no_alpha", System.bp3),
     ("bthp_ones", System.bthp)],
)
def test_determine_input_kind(fixture_name, kind, request):
    assert determine_input_kind(request.getfixturevalue(fixture_name)) == kind


def test_determine_input_kind_fail(example_fail):
    with pytest.raises(UnsupportedTransformationError):
        determine_input_kind(example_fail)


def test_check_all_paths_resolve(request):
    for source_system in System:
        for target_system in System:
            # source_system = System.npol
            # target_system = System.bpb
            if source_system != target_system:  # don't need to convert if the input and output systems are the same
                try:
                    path = get_transform_path(source_system, target_system)
                except UnsupportedTransformationError:
                    pass  # skip any time the transform shouldn't be possible
                else:
                    print(source_system, target_system, path)
                    result = resolve(request.getfixturevalue(f"{source_system}_ones"),
                                     target_system,
                                     out_angles=[0]*u.deg)
                    assert isinstance(result, NDCollection)


def test_btbr_to_npol_missing_out_angles(btbr_ones):
    with pytest.raises(ValueError):
        resolve(btbr_ones, "npol")


def test_imax_effect(mzp_data):
    result = resolve(mzp_data, "mzp", imax_effect=True)
    assert isinstance(result, NDCollection)
    for key in result.keys():
        assert np.sum(result[key].data * mzp_data[key].data) != 0


def test_imax_effect_unsupported_transformation_input(bpb_data):
    with pytest.raises(UnsupportedTransformationError):
        result = resolve(bpb_data, "MZP", imax_effect=True)
        assert isinstance(result, NDCollection)


def test_mzp_to_mzp_is_constant(mzp_ones_other_order):
    result = resolve(mzp_ones_other_order, "mzp", imax_effect=False)
    assert result == mzp_ones_other_order


def test_mzp_to_npol_as_mzp_is_constant(mzp_ones):
    result = resolve(mzp_ones, "npol", out_angles=[-60, 0, 60]*u.degree, imax_effect=False)
    assert np.allclose(result['60.0 deg'].data, mzp_ones["P"].data)
    assert np.allclose(result['-60.0 deg'].data, mzp_ones["M"].data)
    assert np.allclose(result['0.0 deg'].data, mzp_ones["Z"].data)
