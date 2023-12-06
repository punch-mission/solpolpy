import os
from ndcube import NDCollection
import pytest

from solpolpy.instruments import load_data, get_instrument_mask
from solpolpy.errors import TooFewFilesError


def test_load_data():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)
    assert isinstance(out, NDCollection)


def test_load_data_fails_with_one_file():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts"]
    with pytest.raises(TooFewFilesError):
        out = load_data(file_list)


def test_get_instrument_mask():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)

    mask = get_instrument_mask(out['angle_1'].meta)

    assert mask.dtype == bool
    assert mask.shape == (2048, 2048)


def test_use_instrument_mask_overrides_in_load():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list, use_instrument_mask=True)

    assert out['angle_1'].mask.dtype == bool
    assert out['angle_1'].mask.shape == (2048, 2048)

