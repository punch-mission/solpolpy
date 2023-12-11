import os

from solpolpy.instruments import load_data
from solpolpy.plotting import get_colormap_str, plot_collection


def test_get_colormap_str_stereo():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)
    found = get_colormap_str(out['angle_1'].meta)

    assert found == "stereocor2"


def test_get_colormap_str_lasco():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"lasco_0.fts",
               path_to_test_files+"lasco_-60.fts",
               path_to_test_files+"lasco_+60.fts"]
    out = load_data(file_list)
    found = get_colormap_str(out['angle_1'].meta)

    assert found == "soholasco2"


def test_get_colormap_str_unknown():
    meta = {'INSTRUME': 'unknown'}
    found = get_colormap_str(meta)
    assert found == "soholasco2"


def test_plot_collection_runs():
    """ it's hard to test matplotlib things... so this at least makes sure it runs without error"""
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)

    plot_collection(out, show_colorbar=True)
