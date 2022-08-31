"""
Created by Bryce M. Walbridge on 6/02/2022.

The module explored definitions of math equation illustrated by the article,
Three-polarizer Treatment of Linear Polarization in Coronagraphs and
Heliospheric Imagers by Craig E. DeForest, Daniel B. Seaton, and
Matthew J. West. The overall goal is to convert values from three polarized
images (B, pB pB') or Stokes (I, Q, U) representations of linear
polarization from polarizer triplet data. This based on the formulae proposed
in this article.
"""

import numpy as np


def del_Bi(delta_c_B, B_0, Bi):
    """The module calculates the delta of total brightness for ith angles by using parameters
    of signal independent noise `delta_c_B`, an instrument-specific
constant of proportionality, `B_0`, and using total brightness at the ith angle, `Bi`. Note that,
sqrt(B_0*Bi) in this equation is thePoisson noise term that arises from counting statistics.

    This function takes in a var of `delta_c_B`, `B_0`, `Bi` .

    Parameters
    ----------
    delta_c_B : np.array
    B_0: int/ np.array depending on number of intruments
    Bi : dict

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_Bi`.

    Notes
    ------
    Equation 19 of DeForest et al. 2022
    """
    delta_Bi = np.sqrt(((delta_c_B) ** 2) + (B_0 * Bi))

    return delta_Bi


def del_B(delta_c_B, B_0, B):
    """The module calculates the delta of total brightness using parameters
    of signal independent noise `delta_c_B`, an instrument-specific
constant of proportionality, `B_0`, and using total brightness, `B`.

    This function takes in a var of `delta_c_B`, `B_0`, `B` .

    Parameters
    ----------
    delta_c_B : np.array
    B_0: int/ np.array depending on number of intruments
    B : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_B`.

    Notes
    ------
    Equation 22 of DeForest et al. 2022
    """

    delta_B = (2 / np.sqrt(3)) * np.sqrt((delta_c_B) ** 2 + ((B_0 * B) / 2))

    return delta_B


def pB_error(delta_c_B, B_0, B, theta_i, pB, pB_prime, alpha):
    """The module calculates the delta/error of polarized brightness for ith angles by using parameters
    of signal independent noise `delta_c_B`, an instrument-specific
constant of proportionality, `B_0`, total brightness, `B`, polarized angle, `theta_i`, original polarizer brightnesses,
`pB`, `pB_prime`, and solar position angle of an image point, `alpha`.

    This function takes in a var of `delta_c_B`, `B_0`, `B`, `theta_i`, `pB`, `pB_prime`, and `alpha`.

    Parameters
    ----------
    delta_c_B : np.array
    B_0: int/ np.array depending on number of intruments
    Bi : dict
    theta_i : np.array
    pB : np.array
    pB_prime : np.array
    alpha : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_pB`.

    Notes
    ------
    Equation 26 of DeForest et al. 2022
    """

    part_1 = ((3 / 2) * (delta_c_B) ** 2)
    part_2 = ((3 / 2) * ((B_0 * B) / 2))
    part_3 = (B_0 / 2)

    pB_sum_1 = np.sum(np.cos(2 * (theta_i - alpha)) ** 3)
    pB_sum_2 = np.sum(np.cos(2 * (theta_i - alpha)) ** 2)

    delta_pB = (4 / 3) * np.sqrt(part_1 + part_2
                                 - part_3 * ((pB * (pB_sum_1))
                                             + (pB_prime * (pB_sum_2))))
    return delta_pB


def del_B_through_del_theta(pB, pB_prime, delta_theta):
    """The module calculates the delta of total brightness using parameters
    original polarizer brightnesses, `pB`, `pB_prime`, and the delta of polarized angle, `delta-theta`.

    This function takes in a var of `pB`, `pB_prime`, `delta_theta` .

    Parameters
    ----------
    pB : np.array
    pB_prime: np.array
    delta_theta : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_B`.

    Notes
    ------
    Equation 32 of DeForest et al. 2022
    """

    delta_B = np.sqrt(2 / 3) * np.sqrt((pB ** 2) + (pB_prime ** 2)) * (delta_theta)

    return delta_B


def del_pB(pB_prime, delta_theta):
    """The module calculates the delta of polarized brightness using parameters
    original polarizer brightness, `pB_prime`, and the delta of polarized angle, `delta-theta`.

    This function takes in a var of `pB_prime`, `delta_theta`.

    Parameters
    ----------
    pB_prime: np.array
    delta_theta : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_pB`.

    Notes
    ------
    Equation 37 of DeForest et al. 2022
    """

    delta_pB = 2 * (np.sqrt(2 / 3)) * pB_prime * delta_theta

    return delta_pB


def del_pB_prime(pB, delta_theta):
    """The module calculates the delta of polarized brightness using parameters
    original polarizer brightness, `pB`, and the delta of polarized angle, `delta-theta`.

    This function takes in a var of `pB`, `delta_theta`.

    Parameters
    ----------
    pB: np.array
    delta_theta : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `delta_pB_prime`.

    Notes
    ------
    Equation 38 of DeForest et al. 2022
    """

    delta_pB_prime = 2 * (np.sqrt(2 / 3)) * pB * delta_theta

    return delta_pB_prime


def B_theta_bar_through_B_theta(B_theta, epsilon, B):
    """The module calculates observed polarizer brightness using parameters
    radiance through a polarizer at angle theta, `B_theta`, a small leakage coefficient,
    `epsilon`, and total polarizer brightness, `B`.

    This function takes in a var of `B_theta`, `epsilon`, and `B`.

    Parameters
    ----------
    B_theta: np.array
    epsilon : float
    B : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `B_theta_bar`.

    Notes
    ------
    Equation 39 of DeForest et al. 2022
    """

    B_theta_bar = B_theta + epsilon * B

    return B_theta_bar


def pB_through_epsilon_sum(Bi_bar, B, epsilon, alpha):
    """The module calculates polarized brightness using parameters observed polarizer brightness at an ith angle,
    Bi_bar`, total polarizer brightness, `B`, a small leakage coefficient,
    `epsilon`, and solar position angle of an image point, `alpha`.

    This function takes in a var of `Bi_bar`, `epsiloni`.

    Parameters
    ----------
    Bi_bar : dict
    B : np.array
    epsilon : float
    alpha : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `pB`.

    Notes
    ------
    Equation 40 of DeForest et al. 2022
    """

    pB = (-4 / 3) * np.sum([((Bi_bar - epsilon * B) * ith_polarizer_brightness)
                            * np.cos(2 * (ith_angle - alpha))
                            for ith_angle,
                                ith_polarizer_brightness in Bi_bar.items()], axis=0)
    return pB


def B_through_epsilon_sum(Bi_bar, epsiloni):
    """The module calculates total polarizer brightness using parameters observed polarizer brightness at an ith angle,
    Bi_bar`, a small leakage coefficient summed at different observed polarizer brightnesses, `epsiloni`.

    This function takes in a var of `Bi_bar`, `epsiloni`.

    Parameters
    ----------
    Bi_bar: dict
    epsiloni : np.array

    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `B`.

    Notes
    ------
    Equation 41 of DeForest et al. 2022
    """

    numerator = np.sum(Bi_bar)
    denom = 1 + ((2 / 3) * np.sum(epsiloni))

    B = (2 / 3) * (numerator / denom)

    return B
