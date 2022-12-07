import numpy as np
import astropy.units as u
import pytest
from pytest import fixture

import solpolpy.polarizers as pol
#
# #Calculate M,Z,P using B_theta equation substituting theta for theta_MZP.
# @pytest.mark.parametrize("BR, BT, theta_MZP, alpha, expected",
#                         [(0, 1, -np.pi/3, np.pi/3, 0.75),
#                          (0, 1, 0, np.pi/3, 0.75),
#                          (0, 1, np.pi/3, np.pi/3, 0)])
# def test_answer_9(BR, BT, theta_MZP, alpha, expected):
#     assert np.allclose(eqns.B_theta(BR, BT, theta_MZP, alpha), expected)


@fixture
def mzp_zeros():
    return {np.pi/3 * u.radian: np.array([0]),
            0*u.radian: np.array([0]),
            -np.pi/3 * u.radian: np.array([0]),
            "alpha": np.array([0])*u.radian}


def test_mzp_to_bpb_zeros(mzp_zeros):
    actual = pol.mzp_to_bpb(mzp_zeros)
    expected = {"B": np.zeros(1), "pB": np.zeros(1)}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])

@fixture
def mzp_zeros():
    return {np.pi/3 * u.radian: np.array([0]),
            0*u.radian: np.array([0]),
            -np.pi/3 * u.radian: np.array([0])}
def test_any_to_mzp_zeros(mzp_zeros):
    actual = pol.any_to_mzp(mzp_zeros)
    expected = {"Bm": np.zeros(1), "Bz": np.zeros(1), "Bp": np.zeros(1)}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])


@fixture
def mzp_ones():
    return {np.pi/3 * u.radian: np.array([1]),
            0*u.radian: np.array([1]),
            -np.pi/3 * u.radian: np.array([1]),
            "alpha": np.array([0])*u.radian}


def test_mzp_to_bpb_ones(mzp_ones):
    actual = pol.mzp_to_bpb(mzp_ones)
    expected = {"B": np.array([2]), "pB": np.zeros(1), "alpha": mzp_ones['alpha']}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])


@fixture
def mzp_ones():
    return {np.pi / 3 * u.radian: np.array([1]),
            0 * u.radian: np.array([1]),
            -np.pi / 3 * u.radian: np.array([1])}

def test_npol_to_mzp_ones(mzp_ones):
    actual = pol.npol_to_mzp(mzp_ones)
    expected = {"Bm": np.ones(1), "Bz": np.ones(1), "Bp": np.ones(1)}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])

@fixture
def mzp_fourp():
    return {0*u.radian: np.array([1]),
            np.pi/4*u.radian: np.array([1]),
            np.pi/2*u.radian: np.array([1]),
            3*np.pi/4*u.radian: np.array([1]),}

def test_npol_to_mzp_fourp(mzp_fourp):
    actual = pol.npol_to_mzp(mzp_fourp)
    expected = {"Bm": np.array([4/3]), "Bz": np.array([4/3]), "Bp": np.array([4/3])}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])


@fixture
def mzp_no_alpha():
    return {np.pi/3 * u.radian: np.array([0]),
            0*u.radian: np.array([0]),
            -np.pi/3 * u.radian: np.array([0])}

def test_mzp_to_bpb_missing_alpha_fails(mzp_no_alpha):
    """ This conversion without alpha should fail"""
    with pytest.raises(ValueError):
        pol.mzp_to_bpb(mzp_no_alpha)

@fixture
def mzp_ones():
    return {np.pi/3 * u.radian: np.array([1]),
            0*u.radian: np.array([1]),
            -np.pi/3 * u.radian: np.array([1]),
            "alpha": np.array([30])*u.degree}

def test_mzp_to_bp3(mzp_ones):
    actual = pol.mzp_to_bp3(mzp_ones)
    expected = {"B": np.array([2]), "pB": np.array([0]), "pBp": np.array([0])}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])

def test_mzp_to_stokes(mzp_ones):
    actual = pol.mzp_to_stokes(mzp_ones)
    expected = {"Bi": np.array([2]), "Bq": np.array([0]), "Bu": np.array([0])}
    for k, v in expected.items():
        assert np.allclose(actual[k], expected[k])

# #Calculate B_theta using BR, BT
# @pytest.mark.parametrize("BT, BR, alpha, theta, expected",
#                             [(0, 1, np.pi/3, 0, 0.25),
#                             (1, 0, np.pi/3, 0, 0.75)])
# def test_B_theta(BT, BR, alpha, theta, expected):
#     assert np.allclose(eqns.B_theta(BR, BT, theta, alpha) ,expected)
#
# #Calculate first pB, then B, and expected = (pB, B)
# @pytest.mark.parametrize("BT, BR, expected",
#                             [(0, 1, (-1, 1)),
#                             (1, 0, (1, 1))])
# def test_getting_expected_pB_B_values(BT,BR, expected):
#     assert eqns.inverse_of_radial_tangental_radiance(BT,BR) == expected
#
# #Calculate B_theta using P, pB
# @pytest.mark.parametrize("B, pB, alpha, theta, expected",
#                             [(1, -1, np.pi/3, 0.0, 0.25),
#                             (1, 1, np.pi/3, 0.0, 0.75)])
# def test_B_Theta_through_eqs_4(B, pB, alpha, theta, expected):
#     assert np.allclose(eqns.B_theta_double_angle_formula_substitution(
#                         B, pB, theta, alpha), expected)
#
#
# #Calculate M,Z,P when brightness is completely tangential, alpha=60, theta=angle between reference angle and M,Z,P.
#
#
#
#
# #Using equation 9 in the research article to calculate
# #the Brightness.
# @pytest.mark.parametrize("MZP_Brightness, expected",
#                         [({-np.pi/3: eqns.B_theta(0, 1, -np.pi/3, np.pi/3) ,
#                                   0: eqns.B_theta(0, 1, 0, np.pi/3) ,
#                             np.pi/3: eqns.B_theta(0, 1, np.pi/3, np.pi/3)}, 1.0)])
# def test_B_is_expected(MZP_Brightness, expected):
#     assert np.allclose(eqns.B_MZP_sum(MZP_Brightness), expected)
#
# #Using equation 7 in the research article to calculate
# #the polarized Brightness.
# @pytest.mark.parametrize("MZP_Brightness, alpha, expected",
#                         [({-np.pi/3: eqns.B_theta(0, 1, -np.pi/3, np.pi/3) ,
#                                   0: eqns.B_theta(0, 1, 0, np.pi/3) ,
#                             np.pi/3: eqns.B_theta(0, 1, np.pi/3, np.pi/3)}, np.pi/3, 1.0),
#                          (({-np.pi/3: (1/(np.cos(2*(-np.pi/3 -np.pi/3)))) ,
#                            0: (1/(np.cos(2*(0 -np.pi/3)))) ,
#                            np.pi/3: (1/(np.cos(2*(np.pi/3 -np.pi/3))))}),
#                            np.pi/3, -4)])
# def test_pB_is_expected(MZP_Brightness, alpha, expected):
#     print(eqns.pB_MZP_sum(MZP_Brightness, alpha))
#     assert np.allclose(eqns.pB_MZP_sum(MZP_Brightness, alpha), expected)
#
# #TODO replace np.nana with a number
# # #Testing equation 5
# # @pytest.mark.parametrize("B, B_theta, theta, alpha, expected",
# #                         [(0, 1.0, np.pi/3, np.pi/4, 0)])
# # def test_pB(B, B_theta, theta, alpha, expected):
# #     print(eqns.pB(B, B_theta, theta, alpha))
# #     assert np.allclose(eqns.pB(B, B_theta, theta, alpha), expected)
#
# # #TODO replace np.nana with a number
# #testing equation 6 ????????
# @pytest.mark.parametrize("B, Bi, alpha, expected",
#                         [(1.0, ({0 : np.zeros(()),
#                         120 : np.zeros(()),
#                         240 : np.zeros(())}),
#                         np.pi/4, 1.2)])
# def test_pB_through_sum(B, Bi, alpha, expected):
#     print(eqns.pB_through_sum(B, Bi, alpha))
#     assert np.allclose(eqns.pB_through_sum(B, Bi, alpha), expected, atol= 0.1)
#
# # #testing equation 8
# @pytest.mark.parametrize("pB, B_theta, theta, alpha, expected",
#                         [(1.0, 1.0, np.pi/2, np.pi/4, 2.0)])
# def test_polarized_B(pB, B_theta, theta, alpha, expected):
#     print(eqns.polarized_B(pB, B_theta, theta, alpha))
#     assert np.allclose(eqns.polarized_B(pB, B_theta, theta, alpha), expected)
#
# #TODO change zero to another integer to test np.array
# #testing equation 10
# @pytest.mark.parametrize("Bi, alpha, expected",
#                         [(({0 : np.zeros((2,2)),
#                         120 : np.zeros((2,2)),
#                         240 : np.zeros((2,2))}), 0, 0)])
# def test_pB_prime(Bi, alpha, expected):
#     print(eqns.pB_prime(Bi, alpha))
#     assert np.allclose(eqns.pB_prime(Bi, alpha), expected)
#
# #testing equation 11
# @pytest.mark.parametrize("B, pB, pB_prime, theta, alpha, expected",
#                         [(5.0, 1.0, 0, 0, 0, 2.0)])
# def test_B_theta_using_pB_prime(B, pB, pB_prime, theta, alpha, expected):
#     print(eqns.B_theta_using_pB_prime(B, pB, pB_prime, theta, alpha))
#     assert np.allclose(eqns.B_theta_using_pB_prime(B, pB, pB_prime, theta, alpha), expected)
#
# #testing equation 12
# @pytest.mark.parametrize("B_z,B_m,B_p, expected",
#                         [(1.0, 1.0, 1.0, 0)])
# def test_Stokes_Q(B_z,B_m,B_p, expected):
#     print(eqns.Stokes_Q(B_z,B_m,B_p))
#     assert np.allclose(eqns.Stokes_Q(B_z,B_m,B_p), expected)
#
# #testing equation 13
# @pytest.mark.parametrize("B_p,B_m, expected",
#                         [(1.0, 1.0, 0)])
# def test_Stokes_U(B_p,B_m, expected):
#     print(eqns.Stokes_U(B_p,B_m))
#     assert np.allclose(eqns.Stokes_U(B_p,B_m), expected)
#
# #testing equation 15
# @pytest.mark.parametrize("pB, pB_prime, alpha, expected",
#                         [(1.0, 0, np.pi/2, np.pi)])
# def test_theta_maximum(pB, pB_prime, alpha, expected):
#     print(eqns.theta_maximum(pB, pB_prime, alpha))
#     assert np.allclose(eqns.theta_maximum(pB, pB_prime, alpha), expected)
#
# #testing equation 16
# @pytest.mark.parametrize("B, pB, pB_prime, expected",
#                         [(1.0, 0, 0, 0)])
# def test_polarization_fraction(B, pB, pB_prime, expected):
#     print(eqns.polarization_fraction(B, pB, pB_prime))
#     assert np.allclose(eqns.polarization_fraction(B, pB, pB_prime), expected)
#
# #testing alpha function
# @pytest.mark.parametrize("array_length, expected",
#                         [(0, 0)])
# def test_alpha_func(array_length, expected):
#     print(eqns.alpha_func(array_length))
#     assert np.allclose(eqns.alpha_func(array_length), expected)