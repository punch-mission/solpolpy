import os

import numpy as np
import pytest

from solpolpy import resolve
from solpolpy.instruments import load_data
from solpolpy.plotting import generate_rgb_image, get_colormap_str, plot_collection


def test_get_colormap_str_stereo():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)
    found = get_colormap_str(out["0.0 deg"].meta)

    assert found == "stereocor2"


def test_get_colormap_str_lasco():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"lasco_0.fts",
               path_to_test_files+"lasco_-60.fts",
               path_to_test_files+"lasco_+60.fts"]
    out = load_data(file_list)
    found = get_colormap_str(out["0.0 deg"].meta)

    assert found == "soholasco2"


def test_get_colormap_str_unknown():
    meta = {"INSTRUME": "unknown"}
    found = get_colormap_str(meta)
    assert found == "soholasco2"


def test_plot_collection_runs():
    """It's hard to test matplotlib things... so this at least makes sure it runs without error"""
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)

    plot_collection(out, show_colorbar=True)


@pytest.fixture
def sample_collection():
    """Creates a sample NDCollection for testing."""
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+"/test_support_files/"
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    loaded_data = load_data(file_list)
    input_data = resolve(loaded_data, "mzpsolar")

    return input_data


def test_generate_rgb_image(sample_collection):
    color_image = generate_rgb_image(sample_collection)

    assert isinstance(color_image, np.ndarray)
    assert color_image.shape[0] == 3
    assert color_image.dtype == np.uint8
    assert np.allclose(np.max(color_image), 255)
