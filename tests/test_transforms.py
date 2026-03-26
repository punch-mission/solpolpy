import astropy.units as u
import astropy.wcs
import numpy as np
import pytest
from astropy.io import fits
from ndcube import NDCollection, NDCube
from pytest import fixture

import solpolpy.transforms as transforms
from solpolpy.errors import InvalidDataError, MissingAlphaError, SolpolpyError
from solpolpy.physics.polarization import MZP_ANGLES
from tests.fixtures import *

wcs = astropy.wcs.WCS(naxis=3)
wcs.ctype = "WAVE", "HPLT-TAN", "HPLN-TAN"
wcs.cunit = "Angstrom", "deg", "deg"
wcs.cdelt = 0.2, 0.5, 0.4
wcs.crpix = 0, 2, 2
wcs.crval = 10, 0.5, 1
wcs.cname = "wavelength", "HPC lat", "HPC lon"


def test_bpb_mzp_zeros(bpb_zeros):
    actual = transforms.bpb_to_mzpsolar(bpb_zeros)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_mzp_ones(bpb_ones):
    actual = transforms.bpb_to_mzpsolar(bpb_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([3 / 4]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([3 / 4]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_bpb_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bpb_btbr_ones(btbr_bpb_ones):
    actual = transforms.bpb_to_btbr(btbr_bpb_ones)
    expected_data = []
    expected_data.append(("Bt", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Br", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_ones():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bt"})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Br"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_bpb_ones(btbr_ones):
    actual = transforms.btbr_to_bpb(btbr_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("pB", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_btbr_mzp_ways(btbr_ones):
    actual_mzp_direct = transforms.btbr_to_mzpsolar(btbr_ones)
    actual_bpb = transforms.btbr_to_bpb(btbr_ones)
    actual_mzp_indirect = transforms.bpb_to_mzpsolar(actual_bpb)
    for k in list(actual_mzp_direct):
        assert np.allclose(actual_mzp_direct[str(k)].data, actual_mzp_indirect[str(k)].data)


def test_mzp_stokes_ones(mzpsolar_ones):
    actual = transforms.mzpsolar_to_stokes(mzpsolar_ones)
    expected_data = []
    expected_data.append(("I", NDCube(np.arange(0, 10, 2)[:, None] * np.ones((1, 5)), wcs=wcs)))
    expected_data.append(("Q", NDCube(np.zeros((5,5)), wcs=wcs)))
    expected_data.append(("U", NDCube(np.zeros((5,5)), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def stokes_ones():
    data_out = []
    data_out.append(("I", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bi"})))
    data_out.append(("Q", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bq"})))
    data_out.append(("U", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bu"})))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_stokes_mzp_ones(stokes_ones):
    actual = transforms.stokes_to_mzpsolar(stokes_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([(-np.sqrt(3) + 1) / 4]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([(np.sqrt(3) + 1) / 4]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_stokes_to_mzpsolar_returns_quantity_alpha_plane(stokes_ones):
    actual = transforms.stokes_to_mzpsolar(stokes_ones)
    np.testing.assert_allclose(actual["alpha"].data, np.full_like(actual["M"].data, np.pi / 2), rtol=1e-12, atol=1e-12)


def test_mzp_bp3_missing_alpha_errors(mzpsolar_ones):
    with pytest.raises(MissingAlphaError):
        transforms.mzpsolar_to_bp3(mzpsolar_ones)


@fixture()
def bp3_ones():
    data_out = []
    data_out.append(("B", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "B"})))
    data_out.append(("pB", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pB"})))
    data_out.append(("pBp", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "pBp"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_bp3_mzp_ones(bp3_ones):
    actual = transforms.bp3_to_mzpsolar(bp3_ones)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([(3 + np.sqrt(3)) / 4]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([0.0]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([(3 - np.sqrt(3)) / 4]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


@fixture()
def btbr_ones_mzp():
    data_out = []
    data_out.append(("Bt", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Bt"})))
    data_out.append(("Br", NDCube(np.array([1]), wcs=wcs, meta={"POLAR": "Br"})))
    data_out.append(("alpha", NDCube(np.array([0]) * u.degree, wcs=wcs)))
    return NDCollection(data_out, meta={}, aligned_axes="all")


def test_btbr_mzp_ones(btbr_ones_mzp):
    actual = transforms.btbr_to_mzpsolar(btbr_ones_mzp)
    expected_data = []
    expected_data.append(("M", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("Z", NDCube(np.array([1]), wcs=wcs)))
    expected_data.append(("P", NDCube(np.array([1]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_bp3_bthp_ones(bp3_ones):
    actual = transforms.bp3_to_bthp(bp3_ones)
    expected_data = []
    expected_data.append(("B", NDCube(np.array([1]), wcs=wcs)))
    expected_theta = transforms.wrap_linear_polarization(5 * np.pi / 8 * u.radian).to_value(u.radian)
    expected_data.append(("theta", NDCube(np.array([expected_theta]), wcs=wcs)))
    expected_data.append(("p", NDCube(np.array([np.sqrt(2)]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_btbr_npol_ones(btbr_ones):
    actual = transforms.btbr_to_npol(btbr_ones, [0, 120, 240] * u.degree)
    expected_data = [(str(0 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     (str(120 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     (str(240 * u.degree), NDCube(np.array([1]), wcs=wcs)),
                     ("alpha", NDCube(np.array([0]) * u.radian, wcs=wcs))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)


def test_mzp_to_npol_custom():
    """M, Z, P = 0, 1, 0 conversion"""
    input_data = NDCollection(
        [("P", NDCube(np.array([[0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
         ("Z", NDCube(np.array([[1]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
         ("M", NDCube(np.array([[0]]), wcs=wcs, meta={"POLAR": -60 * u.degree}))],
        meta={}, aligned_axes="all")
    actual = transforms.mzpsolar_to_npol(input_data, out_angles=[0, 45, 90] * u.degree)
    expected_data = [(str(0 * u.degree), NDCube(np.array([1]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                     (str(45 * u.degree), NDCube(np.array([1 / 3]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
                     (str(90 * u.degree), NDCube(np.array([-1 / 3]), wcs=wcs, meta={"POLAR": 90 * u.degree}))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)

def test_mzp_to_npol_many_angles():
    """M, Z, P = 0, 1, 0 conversion"""
    input_data = NDCollection(
        [("P", NDCube(np.array([[1, 2], [3, 4]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
         ("Z", NDCube(np.array([[5, 6], [7, 8]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
         ("M", NDCube(np.array([[9, 10], [11, 12]]), wcs=wcs, meta={"POLAR": -60 * u.degree}))],
        meta={}, aligned_axes="all")
    actual = transforms.mzpsolar_to_npol(input_data, out_angles= np.stack([[0, 45], [-60, -15], [60, 105]]) * u.degree)

    expected_keys = [str(22.0), str(-38.0), str(82.0)]

    for k in range(3):
        assert list(actual)[k] == expected_keys[k]

def test_npol_to_mzp_many_angles():
    input_data = NDCollection(
        [("45", NDCube(np.array([[1, 2, 13], [3, 4, 14]]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
         ("20", NDCube(np.array([[5, 6, 15], [7, 8, 16]]), wcs=wcs, meta={"POLAR": 20 * u.degree})),
         ("55", NDCube(np.array([[9, 10, 17], [11, 12, 18]]), wcs=wcs, meta={"POLAR": 55 * u.degree}))],
        meta={}, aligned_axes="all")

    phi = np.stack([np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]),
                np.array([[15.0, 25.0, 35.0], [45.0, 55.0, 65.0]]),
                np.array([[5.0, 15.0, 25.0], [35.0, 45.0, 55.0]])]) * u.degree

    actual = transforms.npol_to_mzpsolar(input_data, in_angles= phi)

    expected_keys = ["M", "Z", "P"]

    mzp_angles = np.array([-60, 0, 60]) * u.degree

    conv_matrix = np.array([[(4 * np.cos(phi[j] - mzp_angles[i]) ** 2 - 1)
          for i in range(3)] for j in range(3)]) / 3

    expected = np.matmul(np.linalg.inv(np.moveaxis(conv_matrix, (0, 1), (-2, -1))),
        np.stack([input_data[k].data for k in input_data], axis=-1)[..., None])

    np.testing.assert_allclose(np.stack([actual[k].data for k in ["M", "Z", "P"]], axis=-1),
        expected[..., 0], rtol=1e-12)

    for k in range(3):
        assert list(actual)[k] == expected_keys[k]

def test_npol_to_mzp_phi_1D():
    input_data = NDCollection(
        [("45", NDCube(np.array([[1, 2, 13], [3, 4, 14]]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
         ("20", NDCube(np.array([[5, 6, 15], [7, 8, 16]]), wcs=wcs, meta={"POLAR": 20 * u.degree})),
         ("55", NDCube(np.array([[9, 10, 17], [11, 12, 18]]), wcs=wcs, meta={"POLAR": 55 * u.degree}))],
        meta={}, aligned_axes="all")

    phi = np.array([10.0, 25.0, 40.0]) * u.degree
    mzp_angles = np.array([-60, 0, 60]) * u.deg

    actual = transforms.npol_to_mzpsolar(input_data, in_angles= phi)

    conv_matrix = ((4 * np.cos(phi[:, None] - mzp_angles[None, :]) ** 2 - 1) / 3)

    expected = np.matmul(np.linalg.inv(conv_matrix), np.stack([input_data[k].data
                                for k in input_data], axis=-1)[..., None])

    np.testing.assert_allclose(np.stack([actual[k].data for k in ["M", "Z", "P"]],
                                        axis=-1), expected[..., 0])


def test_direct_transform_requires_out_angles(btbr_ones):
    with pytest.raises(InvalidDataError, match="Out angles"):
        transforms.btbr_to_npol(btbr_ones)


def test_direct_transform_infers_in_angles_from_metadata(npol_ones):
    actual = transforms.npol_to_mzpsolar(npol_ones)
    assert list(actual.keys())[:3] == ["M", "Z", "P"]


def test_npol_to_mzpsolar_metadata_inference_uses_solar_angles_only():
    input_data = NDCollection(
        [
            ("45", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 45 * u.degree, "POLAROFF": 5 * u.degree})),
            ("20", NDCube(np.array([[2.0]]), wcs=wcs, meta={"POLAR": 20 * u.degree, "POLAROFF": 5 * u.degree})),
            ("55", NDCube(np.array([[3.0]]), wcs=wcs, meta={"POLAR": 55 * u.degree, "POLAROFF": 5 * u.degree})),
        ],
        meta={},
        aligned_axes="all",
    )

    inferred = transforms.npol_to_mzpsolar(input_data)
    explicit = transforms.npol_to_mzpsolar(input_data, in_angles=np.array([45.0, 20.0, 55.0]) * u.degree)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(inferred[key].data, explicit[key].data, rtol=1e-12, atol=1e-12)


def test_npol_to_mzpsolar_ignores_polaroff_for_solar_referenced_angles():
    input_data = NDCollection(
        [
            ("45", NDCube(np.array([[1.0]]), wcs=wcs, meta={"POLAR": 45 * u.degree, "POLAROFF": 90 * u.degree})),
            ("20", NDCube(np.array([[2.0]]), wcs=wcs, meta={"POLAR": 20 * u.degree, "POLAROFF": 90 * u.degree})),
            ("55", NDCube(np.array([[3.0]]), wcs=wcs, meta={"POLAR": 55 * u.degree, "POLAROFF": 90 * u.degree})),
        ],
        meta={},
        aligned_axes="all",
    )

    inferred = transforms.npol_to_mzpsolar(input_data)
    explicit = transforms.npol_to_mzpsolar(input_data, in_angles=np.array([45.0, 20.0, 55.0]) * u.degree)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(inferred[key].data, explicit[key].data, rtol=1e-12, atol=1e-12)


def _roundtrip_collection(values_by_key, alpha=None):
    cubes = []
    for key, value in values_by_key.items():
        cubes.append((key, NDCube(np.asarray(value, dtype=float), wcs=wcs, meta={"POLAR": key})))
    if alpha is not None:
        cubes.append(("alpha", NDCube(np.asarray(alpha) * u.radian, wcs=wcs, meta={"POLAR": "alpha"})))
    return NDCollection(cubes, meta={}, aligned_axes="all")


def test_mzpsolar_bp3_roundtrip_matches_input():
    alpha = np.array([[0.1, 0.2], [0.3, 0.4]])
    mzp = NDCollection(
        [
            ("M", NDCube(np.array([[2.0, 1.5], [1.0, 0.5]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
            ("Z", NDCube(np.array([[0.5, 1.0], [1.5, 2.0]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
            ("P", NDCube(np.array([[1.25, 0.75], [1.75, 2.25]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
            ("alpha", NDCube(alpha * u.radian, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )

    bp3 = transforms.mzpsolar_to_bp3(mzp)
    roundtrip = transforms.bp3_to_mzpsolar(bp3)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(roundtrip[key].data, mzp[key].data, rtol=1e-12, atol=1e-12)


def test_bp3_to_mzpsolar_matches_explicit_equation():
    B = np.array([[2.0, 3.0], [4.0, 5.0]])
    pB = np.array([[0.3, -0.1], [0.2, -0.4]])
    pBp = np.array([[0.5, 0.2], [-0.3, 0.1]])
    alpha = np.array([[0.0, 0.1], [0.2, -0.1]])
    bp3 = _roundtrip_collection({"B": B, "pB": pB, "pBp": pBp}, alpha=alpha)

    actual = transforms.bp3_to_mzpsolar(bp3)

    for key, angle in zip(["M", "Z", "P"], MZP_ANGLES, strict=False):
        expected = 0.5 * (
            B
            - pB * np.cos(2 * (angle.to_value(u.radian) - alpha))
            - pBp * np.sin(2 * (angle.to_value(u.radian) - alpha))
        )
        np.testing.assert_allclose(actual[key].data, expected, rtol=1e-12, atol=1e-12)


def test_mzpsolar_to_bpb_matches_equations_7_and_9():
    alpha = np.array([[0.0, 0.1], [0.2, -0.1]])
    mzp = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
            ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
            ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
            ("alpha", NDCube(alpha * u.radian, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )

    actual = transforms.mzpsolar_to_bpb(mzp)

    polarizer_stack = np.stack([mzp[key].data for key in ["M", "Z", "P"]], axis=0)
    expected_B = (2.0 / 3.0) * np.sum(polarizer_stack, axis=0)
    expected_pB = (-4.0 / 3.0) * np.sum(
        [
            data * np.cos(2 * (angle.to_value(u.radian) - alpha))
            for data, angle in zip(polarizer_stack, MZP_ANGLES, strict=False)
        ],
        axis=0,
    )

    np.testing.assert_allclose(actual["B"].data, expected_B, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(actual["pB"].data, expected_pB, rtol=1e-12, atol=1e-12)


def test_btbr_to_mzpsolar_matches_equations_1_and_3():
    alpha = np.array([[0.0, 0.1], [0.2, -0.1]])
    Bt = np.array([[2.0, 3.0], [4.0, 5.0]])
    Br = np.array([[1.5, 2.5], [3.5, 4.5]])
    btbr = NDCollection(
        [
            ("Bt", NDCube(Bt, wcs=wcs, meta={"POLAR": "Bt"})),
            ("Br", NDCube(Br, wcs=wcs, meta={"POLAR": "Br"})),
            ("alpha", NDCube(alpha * u.radian, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )

    actual = transforms.btbr_to_mzpsolar(btbr)

    for key, angle in zip(["M", "Z", "P"], MZP_ANGLES, strict=False):
        delta = angle.to_value(u.radian) - alpha
        expected = Bt * np.sin(delta) ** 2 + Br * np.cos(delta) ** 2
        np.testing.assert_allclose(actual[key].data, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    ("forward", "backward", "keys", "alpha_required"),
    [
        (transforms.bpb_to_btbr, transforms.btbr_to_bpb, ["B", "pB"], True),
        (transforms.bpb_to_mzpsolar, transforms.mzpsolar_to_bpb, ["B", "pB"], True),
        (transforms.mzpsolar_to_stokes, transforms.stokes_to_mzpsolar, ["M", "Z", "P"], False),
        (transforms.mzpsolar_to_bp3, transforms.bp3_to_mzpsolar, ["M", "Z", "P"], True),
    ],
)
def test_forward_backward_roundtrip_is_systematic(forward, backward, keys, alpha_required):
    alpha = np.array([[0.05, 0.2], [0.35, 0.5]])

    if keys == ["B", "pB"]:
        original = _roundtrip_collection(
            {"B": np.array([[4.0, 5.0], [6.0, 7.0]]), "pB": np.array([[1.0, 0.5], [0.25, -0.2]])},
            alpha=alpha,
        )
    elif keys == ["M", "Z", "P"]:
        original = NDCollection(
            [
                ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
                ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
                ("alpha", NDCube(alpha * u.radian, wcs=wcs)),
            ]
            if alpha_required
            else [
                ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
                ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
                ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
            ],
            meta={},
            aligned_axes="all",
        )
    else:
        raise AssertionError("Unhandled system under test")

    transformed = forward(original)
    roundtrip = backward(transformed)

    for key in keys:
        np.testing.assert_allclose(roundtrip[key].data, original[key].data, rtol=1e-12, atol=1e-12)


def test_mzpsolar_bpb_mzpsolar_roundtrip_is_exact_when_pbp_is_zero():
    alpha = np.array([[0.1, 0.2], [0.3, 0.4]])
    bpb = _roundtrip_collection(
        {"B": np.array([[4.0, 5.0], [6.0, 7.0]]), "pB": np.array([[1.0, 0.5], [0.25, -0.2]])},
        alpha=alpha,
    )

    mzp = transforms.bpb_to_mzpsolar(bpb)
    recovered = transforms.mzpsolar_to_bpb(mzp)

    for key in ["B", "pB"]:
        np.testing.assert_allclose(recovered[key].data, bpb[key].data, rtol=1e-12, atol=1e-12)


def test_mzpsolar_bpb_warns_when_input_contains_nonzero_pbp():
    alpha = np.array([[0.0, 0.1], [0.2, -0.1]])
    bp3 = _roundtrip_collection(
        {
            "B": np.array([[2.0, 3.0], [4.0, 5.0]]),
            "pB": np.array([[0.3, -0.1], [0.2, -0.4]]),
            "pBp": np.array([[0.5, 0.2], [-0.3, 0.1]]),
        },
        alpha=alpha,
    )
    mzp = transforms.bp3_to_mzpsolar(bp3)

    with pytest.warns(UserWarning, match=r"assumes pBp = 0.*max\(\|pBp/B\|\)="):
        transforms.mzpsolar_to_bpb(mzp)


def test_mzpsolar_btbr_mzpsolar_roundtrip_is_exact_when_pbp_is_zero():
    alpha = np.array([[0.1, 0.2], [0.3, 0.4]])
    bpb = _roundtrip_collection(
        {"B": np.array([[4.0, 5.0], [6.0, 7.0]]), "pB": np.array([[1.0, 0.5], [0.25, -0.2]])},
        alpha=alpha,
    )

    mzp = transforms.bpb_to_mzpsolar(bpb)
    btbr = transforms.bpb_to_btbr(bpb)
    recovered = transforms.btbr_to_mzpsolar(btbr)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(recovered[key].data, mzp[key].data, rtol=1e-12, atol=1e-12)


def test_mzpsolar_to_btbr_warns_when_input_contains_nonzero_pbp():
    alpha = np.array([[0.0, 0.1], [0.2, -0.1]])
    bp3 = _roundtrip_collection(
        {
            "B": np.array([[2.0, 3.0], [4.0, 5.0]]),
            "pB": np.array([[0.3, -0.1], [0.2, -0.4]]),
            "pBp": np.array([[0.5, 0.2], [-0.3, 0.1]]),
        },
        alpha=alpha,
    )
    mzp = transforms.bp3_to_mzpsolar(bp3)

    with pytest.warns(UserWarning, match=r"assumes pBp = 0.*max\(\|pBp/B\|\)="):
        transforms.mzpsolar_to_bpb(mzp)

    btbr = transforms.bpb_to_btbr(transforms.mzpsolar_to_bpb(mzp))
    recovered = transforms.btbr_to_mzpsolar(btbr)

    differences = [
        np.nanmax(np.abs(recovered[key].data - mzp[key].data)) for key in ["M", "Z", "P"]
    ]
    assert max(differences) > 0


def test_npol_mzpsolar_roundtrip_with_pixelwise_angles():
    phi = np.stack(
        [
            np.array([[10.0, 20.0], [30.0, 40.0]]),
            np.array([[35.0, 45.0], [55.0, 65.0]]),
            np.array([[80.0, 90.0], [100.0, 110.0]]),
        ]
    ) * u.degree

    mzp = NDCollection(
        [
            ("M", NDCube(np.array([[0.8, 1.2], [1.6, 2.0]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
            ("Z", NDCube(np.array([[2.1, 1.7], [1.3, 0.9]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
            ("P", NDCube(np.array([[1.4, 1.1], [0.7, 0.3]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
        ],
        meta={},
        aligned_axes="all",
    )

    npol = transforms.mzpsolar_to_npol(mzp, out_angles=phi)
    roundtrip = transforms.npol_to_mzpsolar(npol, in_angles=phi)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(roundtrip[key].data, mzp[key].data, rtol=1e-11, atol=1e-11)


def test_mzpsolar_npol_roundtrip_ignores_polaroff_in_solar_frame():
    phi = np.array([-60.0, 0.0, 60.0]) * u.degree

    mzp = NDCollection(
        [
            ("M", NDCube(np.array([[0.8, 1.2], [1.6, 2.0]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 90 * u.degree, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[2.1, 1.7], [1.3, 0.9]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 90 * u.degree, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.4, 1.1], [0.7, 0.3]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 90 * u.degree, "POLARREF": "Solar"})),
        ],
        meta={},
        aligned_axes="all",
    )

    npol = transforms.mzpsolar_to_npol(mzp, out_angles=phi)
    roundtrip = transforms.npol_to_mzpsolar(npol)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(roundtrip[key].data, mzp[key].data, rtol=1e-12, atol=1e-12)


def test_reference_angle_changes_mzpsolar_to_npol():
    input_data = NDCollection(
        [
            ("P", NDCube(np.array([[1.0, 2.0], [3.0, 4.0]]), wcs=wcs, meta={"POLAR": 60 * u.degree})),
            ("Z", NDCube(np.array([[5.0, 6.0], [7.0, 8.0]]), wcs=wcs, meta={"POLAR": 0 * u.degree})),
            ("M", NDCube(np.array([[9.0, 10.0], [11.0, 12.0]]), wcs=wcs, meta={"POLAR": -60 * u.degree})),
        ],
        meta={},
        aligned_axes="all",
    )

    nominal = transforms.mzpsolar_to_npol(input_data, out_angles=np.array([0.0, 45.0, 90.0]) * u.degree)
    shifted = transforms.mzpsolar_to_npol(
        input_data,
        out_angles=np.array([0.0, 45.0, 90.0]) * u.degree,
        reference_angle=15 * u.degree,
    )

    assert any(
        not np.allclose(nominal[key].data, shifted[key].data, rtol=1e-12, atol=1e-12)
        for key in [str(0 * u.degree), str(45 * u.degree), str(90 * u.degree)]
    )


def test_reference_angle_changes_npol_to_mzpsolar():
    input_data = NDCollection(
        [
            ("45", NDCube(np.array([[1.0, 2.0], [3.0, 4.0]]), wcs=wcs, meta={"POLAR": 45 * u.degree})),
            ("20", NDCube(np.array([[5.0, 6.0], [7.0, 8.0]]), wcs=wcs, meta={"POLAR": 20 * u.degree})),
            ("55", NDCube(np.array([[9.0, 10.0], [11.0, 12.0]]), wcs=wcs, meta={"POLAR": 55 * u.degree})),
        ],
        meta={},
        aligned_axes="all",
    )

    nominal = transforms.npol_to_mzpsolar(input_data, in_angles=np.array([45.0, 20.0, 55.0]) * u.degree)
    shifted = transforms.npol_to_mzpsolar(
        input_data,
        in_angles=np.array([45.0, 20.0, 55.0]) * u.degree,
        reference_angle=15 * u.degree,
    )

    assert any(
        not np.allclose(nominal[key].data, shifted[key].data, rtol=1e-12, atol=1e-12)
        for key in ["M", "Z", "P"]
    )

@fixture()
def fourpol_ones():
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [(str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0})),
                (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45})),
                (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90})),
                (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}))]
    return NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")


def test_fourpol_to_stokes_ones(fourpol_ones):
    actual = transforms.fourpol_to_stokes(fourpol_ones)
    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs)))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs)))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs)))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)

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

wcs_new = astropy.wcs.WCS(hdr)

@fixture()
def mzp_ones_instru():
    data, _ = np.mgrid[0:5, 0:5]
    input_data = NDCollection(
        [("P", NDCube(data, wcs=wcs_new, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
         ("Z", NDCube(data, wcs=wcs_new, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'})),
         ("M", NDCube(data, wcs=wcs_new, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Instrument'}))],
        meta={}, aligned_axes="all")
    return input_data


def test_mzp_mzp_ones_instru(mzp_ones_instru):
    actual = transforms.mzpinstru_to_mzpsolar(mzp_ones_instru)
    outinstru = transforms.mzpsolar_to_mzpinstru(actual)
    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(outinstru[key].data, mzp_ones_instru[key].data, rtol=1e-12, atol=1e-12)
        assert actual[key].meta["POLARREF"] == "Solar"
        assert u.Quantity(actual[key].meta["POLAROFF"]).to_value(u.degree) == 1.0
        assert u.Quantity(outinstru[key].meta["POLAR"]).to_value(u.degree) in (-60.0, 0.0, 60.0)
        assert u.Quantity(outinstru[key].meta["POLAROFF"]).to_value(u.degree) == 1.0

@fixture()
def mzp_ones_solar():
    input_data = NDCollection(
        [("P", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs_new, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
         ("Z", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs_new, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'})),
         ("M", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs_new, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": 'Solar'}))],
        meta={}, aligned_axes="all")
    return input_data

def test_mzp_mzp_ones_solar(mzp_ones_solar):
    actual = transforms.mzpsolar_to_mzpinstru(mzp_ones_solar)
    recovered = transforms.mzpinstru_to_mzpsolar(actual)
    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(recovered[key].data, mzp_ones_solar[key].data, rtol=1e-12, atol=1e-12)
        assert actual[key].meta["POLARREF"] == "Instrument"

def test_mzp_two_ways(mzp_ones_instru):
    outsolar = transforms.mzpinstru_to_mzpsolar(mzp_ones_instru)
    outinstru = transforms.mzpsolar_to_mzpinstru(outsolar)
    assert np.allclose(outinstru["M"].data, mzp_ones_instru["M"].data, rtol=0.1)
    assert np.allclose(outinstru["Z"].data, mzp_ones_instru["Z"].data, rtol=0.1)
    assert np.allclose(outinstru["P"].data, mzp_ones_instru["P"].data, rtol=0.1)


def test_mzpsolar_to_mzpinstru_matches_equation_45_projection():
    solar = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs_new,
                         meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs_new,
                         meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs_new,
                         meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
        ],
        meta={},
        aligned_axes="all",
    )

    actual = transforms.mzpsolar_to_mzpinstru(solar)

    shape = solar["Z"].data.shape
    lats = transforms.compute_lats(solar["Z"].wcs, shape)
    solar_north = {
        key: transforms.solnorth_from_wcs(solar[key].wcs, shape=shape, precomputed_lats=lats)
        for key in ["M", "Z", "P"]
    }

    source_stack = np.stack([solar[key].data for key in ["M", "Z", "P"]], axis=0)
    source_angles = np.array([-60, 0, 60]) * u.degree
    expected = {}
    for key, nominal in zip(["M", "Z", "P"], source_angles, strict=False):
        alpha_j = nominal + 1 * u.degree
        phi = ((alpha_j - solar_north[key] + 90 * u.degree) % (180 * u.degree)) - 90 * u.degree
        coeffs = [(4 * np.cos((phi - src_angle).to(u.radian).value) ** 2 - 1) / 3 for src_angle in source_angles]
        expected[key] = sum(data * coeff for data, coeff in zip(source_stack, coeffs, strict=False))

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(actual[key].data, expected[key], rtol=1e-12, atol=1e-12)


def test_instrument_frame_polarizer_angles_match_solar_north_convention():
    solar = NDCollection(
        [
            ("M", NDCube(np.ones((2, 2)), wcs=wcs_new, meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("Z", NDCube(np.ones((2, 2)), wcs=wcs_new, meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("P", NDCube(np.ones((2, 2)), wcs=wcs_new, meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
        ],
        meta={},
        aligned_axes="all",
    )

    actual = transforms._instrument_frame_polarizer_angles(solar)
    shape = solar["Z"].data.shape
    lats = transforms.compute_lats(solar["Z"].wcs, shape)
    expected = np.stack(
        [
            transforms.wrap_linear_polarization(
                transforms.solnorth_from_wcs(solar[key].wcs, shape=shape, precomputed_lats=lats)
                + angle
                + 1 * u.degree
            )
            for key, angle in zip(["M", "Z", "P"], MZP_ANGLES, strict=False)
        ]
    )

    np.testing.assert_allclose(actual.to_value(u.degree), expected.to_value(u.degree), rtol=1e-12, atol=1e-12)


def test_mzpinstru_mzpsolar_equation_45_roundtrip_with_dummy_values():
    solar = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs_new,
                         meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs_new,
                         meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs_new,
                         meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
        ],
        meta={},
        aligned_axes="all",
    )

    instru = transforms.mzpsolar_to_mzpinstru(solar)
    recovered = transforms.mzpinstru_to_mzpsolar(instru)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(recovered[key].data, solar[key].data, rtol=1e-12, atol=1e-12)


def test_mzpinstru_mzpsolar_roundtrip_with_nontrivial_data_and_offset():
    instru = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.5], [0.7, 1.9]]), wcs=wcs_new,
                         meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
            ("Z", NDCube(np.array([[0.6, 1.2], [1.8, 0.9]]), wcs=wcs_new,
                         meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
            ("P", NDCube(np.array([[1.7, 0.8], [1.0, 2.1]]), wcs=wcs_new,
                         meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
        ],
        meta={},
        aligned_axes="all",
    )

    solar = transforms.mzpinstru_to_mzpsolar(instru)
    recovered = transforms.mzpsolar_to_mzpinstru(solar)

    for key in ["M", "Z", "P"]:
        np.testing.assert_allclose(recovered[key].data, instru[key].data, rtol=1e-12, atol=1e-12)


def test_reference_angle_changes_mzpsolar_to_mzpinstru():
    solar = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.3], [1.5, 1.7]]), wcs=wcs_new,
                         meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[0.9, 1.2], [1.4, 1.8]]), wcs=wcs_new,
                         meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.6, 0.7], [1.1, 2.0]]), wcs=wcs_new,
                         meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Solar"})),
        ],
        meta={},
        aligned_axes="all",
    )

    nominal = transforms.mzpsolar_to_mzpinstru(solar)
    shifted = transforms.mzpsolar_to_mzpinstru(solar, reference_angle=15 * u.degree)

    assert any(
        not np.allclose(nominal[key].data, shifted[key].data, rtol=1e-12, atol=1e-12)
        for key in ["M", "Z", "P"]
    )


def test_reference_angle_changes_mzpinstru_to_mzpsolar():
    instru = NDCollection(
        [
            ("M", NDCube(np.array([[1.1, 1.5], [0.7, 1.9]]), wcs=wcs_new,
                         meta={"POLAR": -60 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
            ("Z", NDCube(np.array([[0.6, 1.2], [1.8, 0.9]]), wcs=wcs_new,
                         meta={"POLAR": 0 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
            ("P", NDCube(np.array([[1.7, 0.8], [1.0, 2.1]]), wcs=wcs_new,
                         meta={"POLAR": 60 * u.degree, "POLAROFF": 1, "POLARREF": "Instrument"})),
        ],
        meta={},
        aligned_axes="all",
    )

    nominal = transforms.mzpinstru_to_mzpsolar(instru)
    shifted = transforms.mzpinstru_to_mzpsolar(instru, reference_angle=15 * u.degree)

    assert any(
        not np.allclose(nominal[key].data, shifted[key].data, rtol=1e-12, atol=1e-12)
        for key in ["M", "Z", "P"]
    )


def test_uniform_polaroff_propagates_through_reduced_systems():
    solar = NDCollection(
        [
            ("M", NDCube(np.array([[1.1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 7 * u.degree, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[0.9]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 7 * u.degree, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.6]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 7 * u.degree, "POLARREF": "Solar"})),
            ("alpha", NDCube(np.array([[0.1]]) * u.radian, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )

    for collection in [
        transforms.mzpsolar_to_bpb(solar),
        transforms.mzpsolar_to_bp3(solar),
        transforms.mzpsolar_to_stokes(solar),
        transforms.mzpsolar_to_npol(solar, out_angles=np.array([-60.0, 0.0, 60.0]) * u.degree),
    ]:
        for key in collection.keys():
            if key == "alpha":
                continue
            assert u.Quantity(collection[key].meta["POLAROFF"]).to_value(u.degree) == 7.0


def test_mixed_polaroff_warns_and_is_not_silently_collapsed():
    solar = NDCollection(
        [
            ("M", NDCube(np.array([[1.1]]), wcs=wcs, meta={"POLAR": -60 * u.degree, "POLAROFF": 1 * u.degree, "POLARREF": "Solar"})),
            ("Z", NDCube(np.array([[0.9]]), wcs=wcs, meta={"POLAR": 0 * u.degree, "POLAROFF": 2 * u.degree, "POLARREF": "Solar"})),
            ("P", NDCube(np.array([[1.6]]), wcs=wcs, meta={"POLAR": 60 * u.degree, "POLAROFF": 3 * u.degree, "POLARREF": "Solar"})),
            ("alpha", NDCube(np.array([[0.1]]) * u.radian, wcs=wcs)),
        ],
        meta={},
        aligned_axes="all",
    )

    with pytest.warns(UserWarning, match="mixed POLAROFF"):
        bpb = transforms.mzpsolar_to_bpb(solar)

    for key in ["B", "pB"]:
        assert "POLAROFF" not in bpb[key].meta


@fixture()
def npol_degenerate():
    data_out = [("60.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("0.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": 60*u.degree})),
                ("-60.0 deg", NDCube(np.array([[1, 2], [2, 4]]), wcs=wcs, meta={"POLAR": -60*u.degree}))]
    return NDCollection(data_out, meta={}, aligned_axes="all")

def test_npol_degenerate(npol_degenerate):
    with pytest.raises(SolpolpyError, match="Conversion matrix is degenerate"):
        transforms.npol_to_mzpsolar(npol_degenerate, in_angles=[60, 60, -60] * u.degree)


def test_mask_propagation_works_when_none_provided(fourpol_ones):
    actual = transforms.fourpol_to_stokes(fourpol_ones)
    expected = NDCollection(
        [("I", NDCube(np.array([2]), wcs=wcs, mask=None)),
         ("Q", NDCube(np.array([0]), wcs=wcs, mask=None)),
         ("U", NDCube(np.array([0]), wcs=wcs, mask=None))],
        meta={},
        aligned_axes="all")
    for key in list(expected):
        assert np.allclose(actual[key].data, expected[key].data)
        assert actual[key].mask is None


@pytest.mark.parametrize(
    "transform_fn,input_collection",
    [
        (transforms.mzpsolar_to_stokes, "mzpsolar_ones"),
        (transforms.mzpsolar_to_bpb, "mzpsolar_ones_alpha"),
        (transforms.mzpsolar_to_bp3, "mzpsolar_ones_alpha"),
        (transforms.bpb_to_mzpsolar, "bpb_ones"),
        (transforms.btbr_to_mzpsolar, "btbr_ones"),
        (transforms.btbr_to_npol, "btbr_ones"),
        (transforms.bpb_to_btbr, "bpb_ones"),
        (transforms.btbr_to_bpb, "btbr_ones"),
        (transforms.stokes_to_mzpsolar, "stokes_ones"),
        (transforms.fourpol_to_stokes, "fourpol_ones"),
    ],
)
def test_reference_angle_is_inert_for_transforms_that_do_not_use_it(transform_fn, input_collection, request):
    collection = request.getfixturevalue(input_collection)
    kwargs = {}
    if transform_fn is transforms.btbr_to_npol:
        kwargs["out_angles"] = np.array([0.0, 45.0, 90.0]) * u.degree

    nominal = transform_fn(collection, **kwargs)
    shifted = transform_fn(collection, reference_angle=15 * u.degree, **kwargs)

    for key in nominal.keys():
        np.testing.assert_allclose(nominal[key].data, shifted[key].data, rtol=1e-12, atol=1e-12)


def test_mask_propagation_works_mixed_normal_case():
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0}, mask=np.zeros(1, dtype=bool))),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert np.equal(actual[str(k)].mask, expected[str(k)].mask)


def test_mask_propagation_works_all_false_normal_case():
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0}, mask=np.zeros(1, dtype=bool))),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = []
    expected_data.append(("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected_data.append(("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))))
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert np.equal(actual[str(k)].mask, expected[str(k)].mask)


def test_mask_propagation_works_when_not_all_specified(fourpol_ones):
    # set up some data with mixed masks
    wcs = astropy.wcs.WCS(naxis=1)
    wcs.ctype = "ONE"
    wcs.cunit = "deg"
    wcs.cdelt = 0.1
    wcs.crpix = 0
    wcs.crval = 0
    wcs.cname = "ONE"

    data_out = [
        (str(0 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 0})),
        (str(45 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 45}, mask=np.ones(1, dtype=bool))),
        (str(90 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 90}, mask=np.zeros(1, dtype=bool))),
        (str(135 * u.degree), NDCube(data=np.array([1]), wcs=wcs, meta={"POLAR": 135}, mask=np.zeros(1, dtype=bool)))]
    fourpol_ones = NDCollection(key_data_pairs=data_out, meta={}, aligned_axes="all")

    actual = transforms.fourpol_to_stokes(fourpol_ones)

    expected_data = [("I", NDCube(np.array([2]), wcs=wcs, mask=np.ones(1, dtype=bool))),
                     ("Q", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool))),
                     ("U", NDCube(np.array([0]), wcs=wcs, mask=np.ones(1, dtype=bool)))]
    expected = NDCollection(expected_data, meta={}, aligned_axes="all")
    for k in list(expected):
        assert np.allclose(actual[str(k)].data, expected[str(k)].data)
        assert actual[str(k)].mask is None
