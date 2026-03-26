import os

import pytest
from ndcube import NDCollection

from solpolpy.errors import TooFewFilesError
from solpolpy.instruments import construct_mask, get_instrument_mask, load_data


def test_load_data():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)
    actual_keys = list(out)
    expected_keys = ["0.0 deg", "120.0 deg", "240.0 deg"]
    assert actual_keys == expected_keys
    assert isinstance(out, NDCollection)


def test_load_data_fails_with_one_file():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts"]
    with pytest.raises(TooFewFilesError):
        out = load_data(file_list)
        return out

def test_get_instrument_mask():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)

    mask = get_instrument_mask(out["0.0 deg"].meta)

    assert mask.dtype == bool
    assert mask.shape == (2048, 2048)


def test_use_instrument_mask_overrides_in_load():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list, use_instrument_mask=True)

    assert out["0.0 deg"].mask.dtype == bool
    assert out["0.0 deg"].mask.shape == (2048, 2048)


def test_construct_mask_respects_x_y_centers():
    mask = construct_mask(inner_radius=1, outer_radius=9, center_x=3, center_y=1, shape=(4, 6))

    assert mask[1, 3]  # inner radius masks the center pixel
    assert not mask[1, 1]  # a valid annulus pixel stays unmasked
    assert mask[0, 0]  # far edge lies outside the outer radius
