import os

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import StdDevUncertainty, VarianceUncertainty
from astropy.wcs import WCS
from astropy.wcs.utils import add_stokes_axis_to_wcs
from ndcube import NDCollection, NDCube

from solpolpy import resolve
from solpolpy.core import determine_input_kind
from solpolpy.errors import TooFewFilesError
from solpolpy.instruments import get_instrument_mask, load_data, load_punch


def test_load_data():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR + "/test_support_files/"
    file_list = [
        path_to_test_files + "stereo_0.fts",
        path_to_test_files + "stereo_120.fts",
        path_to_test_files + "stereo_240.fts",
    ]
    out = load_data(file_list)
    actual_keys = list(out)
    expected_keys = ["0.0 deg", "120.0 deg", "240.0 deg"]
    assert actual_keys == expected_keys
    assert isinstance(out, NDCollection)


def test_load_data_fails_with_one_file():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR + "/test_support_files/"
    file_list = [path_to_test_files + "stereo_0.fts"]
    with pytest.raises(TooFewFilesError):
        out = load_data(file_list)
        return out


def test_get_instrument_mask():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR + "/test_support_files/"
    file_list = [
        path_to_test_files + "stereo_0.fts",
        path_to_test_files + "stereo_120.fts",
        path_to_test_files + "stereo_240.fts",
    ]
    out = load_data(file_list)

    mask = get_instrument_mask(out["0.0 deg"].meta)

    assert mask.dtype == bool
    assert mask.shape == (2048, 2048)


def test_use_instrument_mask_overrides_in_load():
    TESTDATA_DIR = os.path.dirname(__file__)
    path_to_test_files = TESTDATA_DIR + "/test_support_files/"
    file_list = [
        path_to_test_files + "stereo_0.fts",
        path_to_test_files + "stereo_120.fts",
        path_to_test_files + "stereo_240.fts",
    ]
    out = load_data(file_list, use_instrument_mask=True)

    assert out["0.0 deg"].mask.dtype == bool
    assert out["0.0 deg"].mask.shape == (2048, 2048)


@pytest.fixture
def punch_mosaic():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN"]
    wcs.wcs.crpix = [2, 2]
    wcs.wcs.crval = [1, 1]
    wcs.wcs.cdelt = [0.1, 0.1]
    meta = {
        "OBSLAYR1": "Polar_M",
        "OBSLAYR2": "Polar_Z",
        "OBSLAYR3": "Polar_P",
        "POLARREF": "Solar",
        "DATE-OBS": "2025-09-11T13:02:29",
    }
    data = np.stack([np.full((3, 4), value) for value in (4.0, 2.0, -1.0)])
    return NDCube(
        data, wcs=add_stokes_axis_to_wcs(wcs, 2), meta=meta, uncertainty=StdDevUncertainty(np.full(data.shape, 0.5))
    )


@pytest.mark.parametrize("as_fits", [False, True])
@pytest.mark.parametrize(
    "level,product,labels",
    [
        ("2", "PT", ("M", "Z", "P")),
        ("3", "PA", ("tB", "pB", "pBp")),
        ("3", "PI", ("M", "Z", "P")),
        ("3", "PT", ("tB", "pB", "pBp")),
        ("3", "PA", ("M", "Z", "P")),
    ],
)
def test_load_punch_products(tmp_path, punch_mosaic, level, product, labels, as_fits):
    header = punch_mosaic.wcs.to_header()
    header.update(punch_mosaic.meta)
    header["TYPECODE"] = product
    header["LEVEL"] = level
    header["OBSCODE"] = "M"
    # Reordered planes must follow labels, even for experiment-specific PAM bases.
    data = punch_mosaic.data[::-1].copy()
    for i, label in enumerate(labels[::-1], 1):
        header[f"OBSLAYR{i}"] = f"Polar_{label}"
    ratio = data / 0.5
    ratio[0, 0, 0] = 0
    if as_fits:
        path = tmp_path / "arbitrary_name.fits"
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.CompImageHDU(data, header=header),
                fits.CompImageHDU(ratio, name="Uncertainty array", quantize_level=0),
            ]
        ).writeto(path)
        out = load_punch(path)
    else:
        sigma = np.full(data.shape, 0.5)
        sigma[0, 0, 0] = np.inf
        cube = NDCube(data, wcs=punch_mosaic.wcs, meta=dict(header), uncertainty=StdDevUncertainty(sigma))
        out = load_punch(cube)
    actual_keys = list(out)
    expected_keys = list(labels)

    assert actual_keys == expected_keys
    assert isinstance(out, NDCollection)

    for i, label in enumerate(labels):
        cube = out[label]
        np.testing.assert_array_equal(cube.data, punch_mosaic.data[i])
        np.testing.assert_allclose(cube.uncertainty.array[1:], 0.5)
        assert cube.meta["DATE-OBS"] == punch_mosaic.meta["DATE-OBS"]
        assert cube.meta["LEVEL"] == level
        assert cube.meta["TYPECODE"] == product
        assert cube.wcs.pixel_n_dim == 2
    assert out[labels[-1]].mask[0, 0]
    assert np.isinf(out[labels[-1]].uncertainty.array[0, 0])
    assert determine_input_kind(out) == ("mzpsolar" if labels[0] == "M" else "bp3")


def test_load_punch_mosaic_and_resolve(punch_mosaic):
    mask = np.zeros((3, 4), dtype=bool)
    mask[0, 0] = True
    out = load_punch(punch_mosaic, mask=mask)
    assert tuple(out) == ("M", "Z", "P")
    for cube in out.values():
        assert cube.mask[0, 0]
        np.testing.assert_allclose(cube.wcs.pixel_to_world_values(1, 1), (1, 1))
    result = resolve(out, "bp3")
    np.testing.assert_allclose(result["tB"].data, 2 / 3 * (4 + 2 - 1))
    np.testing.assert_allclose(out["Z"].uncertainty.array, 0.5)
    assert punch_mosaic.meta["OBSLAYR1"] == "Polar_M"
    assert punch_mosaic.mask is None


@pytest.mark.parametrize("as_fits", [False, True])
def test_load_punch_individual_triplet(tmp_path, punch_mosaic, as_fits):
    sources = []
    for i, angle in enumerate((-60, 0, 60)):
        meta = {
            "LEVEL": "1",
            "TYPECODE": ("PM", "PZ", "PP")[i],
            "OBSCODE": "1",
            "POLAR": angle,
            "POLARREF": "Instrument",
            "POLAROFF": 2,
        }
        wcs = punch_mosaic.wcs.dropaxis(2)
        if as_fits:
            header = wcs.to_header()
            header.update(meta)
            path = tmp_path / f"frame{i}.fits"
            fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(punch_mosaic.data[i], header=header)]).writeto(path)
            sources.append(path)
        else:
            sources.append(
                NDCube(punch_mosaic.data[i], wcs=wcs, meta=meta, uncertainty=VarianceUncertainty(np.full((3, 4), 0.25)))
            )
    out = load_punch(sources[::-1])
    assert determine_input_kind(out) == "mzpinstru"
    assert out["M"].meta["POLAR"] == -60 * u.deg
    assert out["M"].meta["POLAROFF"] == 2
    assert [out[key].meta["TYPECODE"] for key in ("M", "Z", "P")] == ["PM", "PZ", "PP"]
    assert all(cube.meta["LEVEL"] == "1" for cube in out.values())
    np.testing.assert_array_equal(out["P"].data, punch_mosaic.data[2])
    if as_fits:
        assert out["M"].uncertainty is None
    else:
        np.testing.assert_allclose(out["M"].uncertainty.array, 0.5)


def test_load_punch_mixed_polarref_rejected(punch_mosaic):
    frames = [
        NDCube(punch_mosaic.data[i], wcs=punch_mosaic.wcs.dropaxis(2), meta={"POLAR": angle, "POLARREF": ref})
        for i, (angle, ref) in enumerate(zip((-60, 0, 60), ("Instrument", "Solar", "Solar")))
    ]

    with pytest.raises(ValueError, match="same POLARREF"):
        load_punch(frames)


def test_load_punch_invalid_pixels_and_source_mask(punch_mosaic):
    punch_mosaic.data[0, 0, 0] = np.nan
    punch_mosaic.uncertainty.array[1, 0, 1] = np.inf
    punch_mosaic.uncertainty.array[1, 1, 0] = -0.5
    punch_mosaic.mask = np.zeros(punch_mosaic.shape, dtype=bool)
    punch_mosaic.mask[2, 0, 2] = True
    out = load_punch(punch_mosaic)
    assert out["M"].mask[0, 0]
    assert out["Z"].mask[0, 1]
    assert out["Z"].mask[1, 0]
    assert out["P"].mask[0, 2]
    assert not out["P"].mask[1, 1]  # Negative signals remain valid.


def test_load_punch_missing_plane_labels_rejected(punch_mosaic):
    del punch_mosaic.meta["OBSLAYR3"]
    with pytest.raises(ValueError, match="OBSLAYR"):
        load_punch(punch_mosaic)


def test_load_punch_duplicate_frames_rejected(punch_mosaic):
    cube = NDCube(punch_mosaic.data[0], wcs=punch_mosaic.wcs.dropaxis(2), meta={"POLAR": -60, "POLARREF": "Solar"})
    with pytest.raises(ValueError, match="Duplicate"):
        load_punch([cube, cube, cube])


@pytest.mark.parametrize("mask", [np.zeros((1, 1), dtype=bool), np.zeros((3, 4))])
def test_load_punch_bad_mask_rejected(punch_mosaic, mask):
    with pytest.raises(ValueError, match="mask"):
        load_punch(punch_mosaic, mask=mask)
