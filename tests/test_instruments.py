import os
from ndcube import NDCollection
from solpolpy.instruments import load_data

def test_load_data():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files=TESTDATA_DIR+'/test_support_files/'
    file_list=[path_to_test_files+"stereo_0.fts",
               path_to_test_files+"stereo_120.fts",
               path_to_test_files+"stereo_240.fts"]
    out = load_data(file_list)
    assert isinstance(out, NDCollection)
