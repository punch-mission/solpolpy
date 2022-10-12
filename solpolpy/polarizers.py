#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic coordinate transforms
"""

import numpy as np


def mzp_to_bpb(input_dict):
    """
    Notes
    ------
    Equation 7 and 9 in Deforest et al. 2022.

    """""
    # TODO: need to check if 3 angles are input.
    # TODO: need to check if separated appropriately if not create quality warning.

    alpha = input_dict['alpha']

    if "alpha" not in input_dict:
        raise ValueError("missing alpha")

    B = (2 / 3) * (np.sum([ith_polarizer_brightness
                           for ith_angle, ith_polarizer_brightness
                           in input_dict.items() if ith_angle != "alpha"], axis=0))

    pB = (-4 / 3) * (np.sum([ith_polarizer_brightness
                             * np.cos(2 * (ith_angle - alpha))
                             for ith_angle, ith_polarizer_brightness
                             in input_dict.items() if ith_angle != "alpha"], axis=0))

    return {"B": B, "pB": pB, "alpha": alpha}


def bpb_to_mzp(input_dict):
    """
    Notes
    ------
    Equation 4 in Deforest et al. 2022.
    """
    alpha = input_dict['alpha']
    B, pB = input_dict['B'], input_dict['pB']
    result = {theta: (1 / 2) * (B - pB * (np.cos(2 * (theta - alpha)))) for theta in np.deg2rad([-60, 0, 60])}
    result['alpha'] = alpha
    return result


def bpb_to_btbr(input_dict):
    """
    Notes
    ------
    Equation 1 and 2 in Deforest et al. 2022.
    """
    B, pB = input_dict['B'], input_dict['pB']
    Br = (B - pB) / 2
    Bt = (B + pB) / 2
    return {'Br': Br, "Bt": Bt}


######################################################################################################################
######################################################################################################################
######################################################################################################################



def inverse_of_radial_tangental_radiance(B_tangential, B_radial):
    """Converts tangential radiance,`B_tangential`, and radial radiance,
     `B_radial`, back into unpolarized brightness,`B`, and Coronal polarized brightness, 
     `pB`.
    
    This function takes in two vars of `B_tangential`, tangential radiance, and 
    `B_radial`,radial radiance. 
    
    Parameters
    ----------
    B_tangential : np.ndarrary
    B_radial : np.ndarrary
    
    Returns
    -------
     float
        The float that is returned is defined to be vars `B`, and `pB`.
    """
    pB = B_tangential - B_radial
    B = B_tangential + B_radial
    return pB, B


def B_theta(B_radial, B_tangential, theta, alpha):
    """Converts Polarizer angle,`theta`, Solar position angle of an image point,
     `alpha`, radial radiance, `B_radial`, and tangential radiance, `B_tangential` 
    into Radiance through a polarizer at angle theta,`B_theta`.
    
    This function takes in four vars of `theta`, `alpha`, `B_radial`,and `B_tangential`
    
    Parameters
    ----------
    theta : np.ndarray
    alpha : np.ndarray
    B_tangential : np.ndarray
    B_radial : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `B_theta`.

    Notes
    ------
    Equation 3 in Deforest et al. 2022.
    """
    return  B_tangential*(np.sin(theta - alpha))**2 + B_radial*(np.cos(theta - alpha))**2


def pB(B, B_theta, theta, alpha):
    """Converts unpolarized brightness,`B`, Radiance through a polarizer at angle theta,`B_theta`,
    Polarizer angle,`theta`, and Solar position angle of an image point, `alpha` into Coronal 
    polarized brightness, `pB`.
    
    This function takes in four vars of `B`, `B_theta`, `theta`,and `alpha`.
    
    Parameters
    ----------
    B : np.ndarray
    B_theta : np.ndarray
    theta : np.ndarray
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.

    Notes
    ------
    Equation 5 in Deforest et al. 2022.

    Sets values to nan when denominator <= 1E-6
    """
    #anywhere where the denom is less than 1e-6 make it a nan
    pB_denom = np.abs(np.cos(2*(theta - alpha)))
    pB = (B-(2*B_theta))/ np.cos(2*(theta - alpha))
    pB[pB_denom < 1e-6] = np.nan
    return pB


def pB_through_sum(B, Bi, alpha):
    """Converts unpolarized brightness,`B`, a dictionary of unpolarized brightnesses,`Bi`,
    Polarizer angle,`theta`, and Solar position angle of an image point, `alpha` into Coronal 
    polarized brightness, `pB`. Uses summation as another way to calculate cornal polarized brightness.
    
    This function takes in four vars of `B`, `B_theta`, `theta`,and `alpha`.
    
    Parameters
    ----------
    B : np.ndarray
    Bi : dict
    theta : np.ndarray
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.

    Notes
    -----
    Implements Equation 6 of Deforest et al. (2022)
    """
    numerator = np.sum([(B - 2 * ith_polarizer_brightness)
                        * np.cos(2*(ith_angle - alpha)) 
                        for ith_angle, 
                        ith_polarizer_brightness in Bi.items()], axis=0)

    denominator = np.sum([np.cos(2*(theta_i-alpha))**2 for theta_i in Bi], axis = 0)

    return numerator / denominator




def polarized_B(pB, B_theta, theta, alpha):
    """Converts Radiance through a polarizer at angle theta,`B_theta`,
    Polarizer angle,`theta`, solar position angle of an image point, `alpha`, and polarized brightness, `pB`, 
    into unpolarized brightness,`B`.
    
    
    This function takes in three vars of `B_theta`, `theta`,and `alpha`.
    
    Parameters
    ----------
    B_theta : np.ndarray
    theta : np.ndarray
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.

    Notes
    ------
    Equation 8 in Deforest et al. 2022.
    """
    return 2*B_theta + (np.cos(2*(theta - alpha)))*pB


def pB_prime(Bi, alpha):
    """Converts total polarizer brightness at an ith angle, `Bi`, and solar position angle of an image point,
    `alpha` into polarized brightness, `pB_prime`
    
    This function takes in four vars of `Bi`, and `alpha`.
    
    Parameters
    ----------
    alpha : np.ndarray
    Bi : dict
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB_prime`.

    Notes
    ------
    Equation 10 in Deforest et al. 2022.
    """
    
    pB_prime = (-4/3)*(np.sum([ith_polarizer_brightness 
                         * np.sin(2*(ith_angle - alpha)) 
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    return pB_prime


def B_theta_using_pB_prime(B, pB, pB_prime, theta, alpha):
    """Converts unpolarized brightness,`B`, original polarized brightnesses,`pB`, `pB_prime`,
    Polarizer angle,`theta`, and solar position angle of an image point, `alpha` into radiance 
    through a polarizer at angle theta, `B_theta`.
    
    This function takes in four vars of `B`, `pB`, `pB_prime`, `theta`, and `alpha`.
    
    Parameters
    ----------
    B : np.ndarray
    pB : np.ndarray
    pB_prime : np.ndarray
    theta : np.ndarray
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `B_theta`.

    Notes
    ------
    Equation 11 in Deforest et al. 2022.
    """

    B_theta = (1/2)*(B - np.cos(2*(theta - alpha))*pB 
                    - np.sin(2*(theta - alpha))*pB_prime)
    
    return B_theta


def Stokes_Q(B_z,B_m,B_p):
    """Converts unpolarized brightness zero,`B_z`, unpolarized brightness plus direction,`B_p`,
     and unpolarized brightness minus direction,`B_m`, into the Stokes parameter, `Q`.
   
    
    This function takes in three vars of `B_z`, `B_m`, and `B_p`.
    
    Parameters
    ----------
    B_z : np.ndarray
    B_m : np.ndarray
    B_p : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `Q`.

    Notes
    ------
    Equation 12 in Deforest et al. 2022.
    """
    return (2/3)*(2*B_z - B_m - B_p)


def Stokes_U(B_p,B_m):
    """Converts unpolarized brightness in the plus direction,`B_p`, unpolarized brightness in the minus direction, `B_p`
    into the Stokes parameter, `U`.
    
    This function takes in two vars of `B_p`, and `B_m`.
    
    Parameters
    ----------
    B_m : np.ndarray
    B_p : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `Q`.

    Notes
    ------
    Equation 13 in Deforest et al. 2022.
    """
    return (2/np.sqrt(3))*(B_p - B_m)


def theta_maximum(pB, pB_prime, alpha):
    """Converts Coronal polarized brightness, `pB`, `pb_prime`,
    and Solar position angle of an image point, `alpha`, into Polarizer angle , `theta_max`. 
    
    This function takes in three vars of `pB`, `pB_prime`, and `alpha`.
    
    Parameters
    ----------
    pB : np.ndarray
    pB_prime : np.ndarray
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `theta_max`.

    Notes
    ------
    Equation 15 in Deforest et al. 2022.
    """
    return ((1/2)*np.arctan(pB_prime/pB)) + np.pi/2 + alpha


def polarization_fraction(B, pB, pB_prime):
    """Converts Coronal polarized brightness, `pB`,`pB_prime`, and
    unpolarized brightness,`B`, into a polarization fraction, `p`. Uses thompson scattering,
    pB = (p)*(B), and previous equation substituion to come to a simpified 
    version for a polarization fraction,  `p`.
    
    This function takes in three vars of `pB`, `pB_prime`, and `B`.
    
    Parameters
    ----------
    pB : np.ndarray
    pB_prime : np.ndarray
    B : np.ndarray

    Returns
    -------
     float
        The float that is returned is defined to be var, `theta_max`.

    Notes
    ------
    Equation 16 in Deforest et al. 2022.
    """
    return np.sqrt((pB**2 + pB_prime**2)) / B
