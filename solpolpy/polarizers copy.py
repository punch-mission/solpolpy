
import numpy as np


def radial_radiance(B, pB):
    """Converts unpolarized brightness,`B`, and Coronal polarized brightness,
     `pB`, into radial radiance.
    
    This function takes in two vars of `B`, Unpolarized brightness, and 
    `pB`,Coronal polarized brightness. 
    
    Parameters
    ----------
    B : np.ndarray 
    pB : np.ndarray
    
    Returns
    -------
    float
        The float that is returned is defined to be var `B_radial`.

    Notes
    ------
    Equation 1 in Deforest et al. 2022.
    """
    return (B - pB) / 2


def tangential_radiance(B, pB):
    """Converts unpolarized brightness,`B`, and Coronal polarized brightness,
     `pB`, into tangential radiance.
    
    This function takes in two vars of `B`, Unpolarized brightness, and 
    `pB`,Coronal polarized brightness. 
    
    Parameters
    ----------
    B : np.ndarrary
    pB : np.ndarrary
    
    Returns
    -------
     float
        The float that is returned is defined to be var `B_tangential`.

    Notes
    ------
    Equation 2 in Deforest et al. 2022.
    """
    return (B + pB)/2


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


def B_theta_double_angle_formula_substitution(B, pB, theta, alpha):
    """Converts unpolarized brightness, `B`, Coronal polarized brightness, `pB`, Polarizer angle,`theta`,
    and solar position angle of an image point, `alpha`, into Radiance through a polarizer at angle theta, `B_theta`. 
    
    This function takes in four vars of `B`, `pB`, `theta`, and `alpha`.
    
    Parameters
    ----------
    B : np.ndarray
    pB : np.ndarray
    theta : np.ndarray
    alpha : np.ndarray

    Returns
    -------
     float
        The float that is returned is defined to be var, `B_theta`.

    Notes
    ------
    Equation 4 in Deforest et al. 2022.
    """
    return (1/2)*(B - pB*(np.cos(2*(theta - alpha))))


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

def pB_MZP_sum(Bi, alpha):
    """Converts Solar position angle of an image point,`alpha`, into Coronal 
    polarized brightness, `pB`. Uses summation as another way to calculate cornal polarized brightness.
    
    This function takes in the var of `alpha`.
    
    Parameters
    ----------
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.
        
    Notes
    ------
    Equation 7 in Deforest et al. 2022.
    """""
    pB = (-4/3)*(np.sum([ith_polarizer_brightness 
                         * np.cos(2*(ith_angle - alpha)) 
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    
    return pB


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


def B_MZP_sum(Bi):
    """Converts unpolarized brightness,`B`, into an average unpolarized brightness.
    Uses summation averaging over the three polarizer positions M,Z, and P.
    
    This function takes in a var of, `Bi`, which is a dictionary of these different polarizer positions.
    
    Parameters
    ----------
    Bi : dict

    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.

    Notes
    ------
    Equation 9 in Deforest et al. 2022.
    """
    
    B = (2/3)*(np.sum([ith_polarizer_brightness  
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    
    return B


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

 
def alpha_func(array_length):
    """Makes an array of alpha values on a mesh grid. This function take in the parameter of
    arrary_length to properly parameterize the created meshgrid.
    
    This function takes in a var of `array_length`.
    
    Parameters
    ----------
    array_length : int
    
    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `alpha`.
    """
    x = np.arange(-array_length//2, array_length//2)
    y = np.arange(-array_length//2, array_length//2)
    xx, yy = np.meshgrid(x, y)
    return np.fliplr(np.rot90(np.arctan2(yy, xx), k=1))
