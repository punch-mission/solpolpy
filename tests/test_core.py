"""Docstring for the test_core.py module.

Created by Bryce M. Walbridge on 6/23/2022.

The module explored definitions of math equation illistrated by the article, 
Three-polarizer Treatment of Linear Polarization in Coronagraphs and 
Heliospheric Imagers by Craig E. DeForester, Daniel B. Seaton, and 
Matthew J. West. The overall goal is to test math definitions to document failures
and successes.

"""

import numpy as np
import pytest
from solpolpy import core as eqns


#Calculate B_theta using BR, BT
@pytest.mark.parametrize("BT, BR, alpha, theta, expected", 
                            [(0, 1, np.pi/3, 0, 0.25), 
                            (1, 0, np.pi/3, 0, 0.75)])
def test_B_theta(BT, BR, alpha, theta, expected): 
    assert np.allclose(eqns.B_theta(BR, BT, theta, alpha) ,expected)

#Calculate first pB, then B, and expected = (pB, B) 
@pytest.mark.parametrize("BT, BR, expected", 
                            [(0, 1, (-1, 1)), 
                            (1, 0, (1, 1))])
def test_getting_expected_pB_B_values(BT,BR, expected):
    assert eqns.inverse_of_radial_tangental_radiance(BT,BR) == expected

#Calculate B_theta using P, pB
@pytest.mark.parametrize("B, pB, alpha, theta, expected", 
                            [(1, -1, np.pi/3, 0.0, 0.25), 
                            (1, 1, np.pi/3, 0.0, 0.75)])
def test_B_Theta_through_eqs_4(B, pB, alpha, theta, expected):
    assert np.allclose(eqns.B_theta_double_angle_formula_substitution(
                        B, pB, theta, alpha), expected)


#Calculate M,Z,P when brightness is completely tangential, alpha=60, theta=angle between reference angle and M,Z,P.
BR = 0
BT = 1
alpha = np.pi/3.
theta_M = -np.pi/3.
theta_Z = 0
theta_P = np.pi/3.
M=eqns.B_theta(BR, BT, 
               theta_M, alpha)
Z=eqns.B_theta(BR, BT, 
               theta_Z, alpha)
P=eqns.B_theta(BR, BT, 
               theta_P, alpha)
MZP_Brightness={theta_M:M , theta_Z:Z , 
                theta_P:P}

#Calculate M,Z,P using B_theta equation substituting theta for theta_MZP.
@pytest.mark.parametrize("BR, BT, theta_MZP, alpha, expected",
                        [(0, 1, -np.pi/3, np.pi/3, 0.75), 
                         (0, 1, 0, np.pi/3, 0.75), 
                         (0, 1, np.pi/3, np.pi/3, 0)])
def test_answer_9(BR, BT, theta_MZP, alpha, expected):
    assert np.allclose(eqns.B_theta(BR, BT, theta_MZP, alpha), expected)



#Using equation 9 in the research article to calculate 
#the Brightness.
@pytest.mark.parametrize("MZP_Brightness, expected",
                        [(MZP_Brightness, 1.0)])
def test_B_is_expected(MZP_Brightness, expected):
    assert np.allclose(eqns.B_MZP_sum(MZP_Brightness), expected)

#Using equation 7 in the research article to calculate 
#the polarized Brightness.
@pytest.mark.parametrize("MZP_Brightness, alpha, expected",
                        [(MZP_Brightness, np.pi/3, 1.0), 
                         (({theta_M: (1/(np.cos(2*(theta_M -alpha)))) , 
                           theta_Z: (1/(np.cos(2*(theta_Z -alpha)))) , 
                           theta_P: (1/(np.cos(2*(theta_P -alpha))))}), 
                           np.pi/3, -4)])
def test_pB_is_expected(MZP_Brightness, alpha, expected):
    print(eqns.pB_MZP_sum(MZP_Brightness, alpha))
    assert np.allclose(eqns.pB_MZP_sum(MZP_Brightness, alpha), expected)
    
