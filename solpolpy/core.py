"""Docstring for the math_equations.py module.

Created by Bryce M. Walbridge on 6/02/2022.

The module explored definitions of math equation illistrated by the article, 
Three-polarizer Treatment of Linear Polarization in Coronagraphs and 
Heliospheric Imagers by Craig E. DeForester, Daniel B. Seaton, and 
Matthew J. West. The overall goal is to convert values from three polarized 
images (B, pB pB') or Stokes (I, Q, U) representations of linear
polarization from polarizer triplet data. This based on the formulae proposed
in this article.

"""
import numpy as np


#Equation 1
def radial_radiance(B, pB):
    '''Converts unpolarized brightness,`B`, and Coronal polarized brightness, 
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
        
    Raises
    ------
    check shapes match.
    '''
    # if B.shape != pB.shape:
    #     raise ValueError(f'B and pB shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'pB.shape = {pB.shape}') 
    
    B_radial = (B - pB)/2
    return B_radial
    
#Equation 2
def tangential_radiance(B, pB):
    '''Converts unpolarized brightness,`B`, and Coronal polarized brightness, 
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
    Raises
    ValueError if shapes do match 
    ------
    
    '''
    # if B.shape != pB.shape:
    #     raise ValueError(f'B and pB shapes do not match. ' 
    #                      f'B.shape = {B.shape} ' 
    #                      f'pB.shape = {pB.shape}')
    
    B_tangential = (B + pB)/2
    return B_tangential


#definiton for inverse function
def inverse_of_radial_tangental_radiance(B_tangential, B_radial):
    '''Converts tangential radiance,`B_tangential`, and radial radiance, 
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
    Raises
    ValueError if shapes do match 
    ------
    
    '''
    # if B_tangential.shape != B_radial.shape:
    #     raise ValueError(f'B_radial and B_tangential shapes do not match. ' 
    #                      f'B_tangential.shape = {B_tangential.shape} '
    #                      f'B_radial.shape = {B_radial.shape}')
    
    pB = B_tangential - B_radial
    B = B_tangential + B_radial
    return pB, B


#definition for electric feild amplitude and Equation 3
def B_theta(B_radial, B_tangential, theta, alpha):
    '''Converts Polarizer angle,`theta`, Solar position angle of an image point, 
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
    Raises
    ValueError if shapes do match 
    ------
    
    '''
    #checks to see if all np.ndarrarys match each other in shape
    # if theta.shape != alpha.shape:
    #     raise ValueError(f'theta and alpha shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif theta.shape != B_radial.shape:
    #     raise ValueError(f'theta and B_radial shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'B_radial.shape = {B_radial.shape}')
        
    # elif theta.shape != B_tangential.shape:
    #     raise ValueError(f'theta and B_tangential shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'B_tangential.shape = {B_tangential.shape}')
        
    # elif alpha.shape != B_radial.shape:
    #     raise ValueError(f'alpha and B_radial shapes do not match. ' 
    #                      f'alpha.shape = {alpha.shape} '
    #                      f'B_radial.shape = {B_radial.shape}')
        
    # elif alpha.shape != B_tangential.shape:
    #     raise ValueError(f'alpha and B_tangential shapes do not match. ' 
    #                      f'alpha.shape = {alpha.shape} '
    #                      f'B_tangential.shape = {B_tangential.shape}')
        
    # elif B_radial.shape != B_tangential.shape:
    #     raise ValueError(f'B_radial and B_tangential shapes do not match. ' 
    #                      f'B_radial.shape = {B_radial.shape} '
    #                      f'B_tangential.shape = {B_tangential.shape}')
        
        
    B_theta = B_tangential*(np.sin(theta - alpha))**2 + B_radial*(np.cos(theta - alpha))**2
    
    return B_theta 


#Equation 4
def B_theta_double_angle_formula_substitution(B, pB, theta, alpha):
    '''Converts unpolarized brightness, `B`, Coronal polarized brightness, `pB`, Polarizer angle,`theta`,
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
    Raises
    ValueError if shapes do match 
    ------
    
    '''
    
    #checks to see if all np.ndarrarys match each other in shape
    # if B.shape != pB.shape:
    #     raise ValueError(f'B and alpha shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'pB.shape = {pB.shape}')
        
    # elif B.shape != theta.shape:
    #     raise ValueError(f'B and theta shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'theta.shape = {theta.shape}')
        
    # elif B.shape != alpha.shape:
    #     raise ValueError(f'B and alpha shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif pB.shape != theta.shape:
    #     raise ValueError(f'pB and theta shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'theta.shape = {theta.shape}')
        
    # elif pB.shape != alpha.shape:
    #     raise ValueError(f'pB and alpha shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif theta.shape != alpha.shape:
    #     raise ValueError(f'theta and alpha shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    B_theta = (1/2)*(B - pB*(np.cos(2*(theta - alpha))))
    
    return B_theta 


#sets a nan for cos being close to zero 
#definition for geting pB with theta and alpha and B_theta. Equation 5
def pB(B, B_theta, theta, alpha):
    '''Converts unpolarized brightness,`B`, Radiance through a polarizer at angle theta,`B_theta`, 
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
    Raises
    ValueError if shapes do match 
    ------
    
    '''
    #checks to see if all np.ndarrarys match each other in shape
    # if B.shape != B_theta.shape:
    #     raise ValueError(f'B and B_theta shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'B_theta.shape = {B_theta.shape}')
        
    # elif B.shape != theta.shape:
    #     raise ValueError(f'B and theta shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'theta.shape = {theta.shape}')
        
    # elif B.shape != alpha.shape:
    #     raise ValueError(f'B and alpha shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif B_theta.shape != theta.shape:
    #     raise ValueError(f'B_theta and theta shapes do not match. ' 
    #                      f'B_theta.shape = {B_theta.shape} '
    #                      f'theta.shape = {theta.shape}')
        
    # elif B_theta.shape != alpha.shape:
    #     raise ValueError(f'B_theta and alpha shapes do not match. ' 
    #                      f'B_theta.shape = {B_theta.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif theta.shape != alpha.shape:
    #     raise ValueError(f'theta and B_tangential shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    #anywhere where the denom is less than 1e-6 make it a nan
    pB_denom = np.abs(np.cos(2*(theta - alpha)))
    pB = (B-(2*B_theta))/ np.cos(2*(theta - alpha))
    pB[pB_denom < 1e-6] = np.nan
    return pB

#Equation 6
def pB_through_sum(B, Bi, alpha):
    '''Converts unpolarized brightness,`B`, a dictionary of unpolarized brightnesses,`Bi`, 
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
    Raises
    ?
    ------
    
    '''
#     alpha = np.zeros((2,2))
#     B = np.zeros((2,2))
#     Bi = {0 : np.zeros((2,2)), 
#         120 : np.zeros((2,2)), 
#         240 : np.zeros((2,2))}

    Numerator = np.sum([(B - 2* ith_polarizer_brightness)  
                        * np.cos(2*(ith_angle - alpha)) 
                        for ith_angle, 
                        ith_polarizer_brightness in Bi.items()], axis=0)

    Denom = np.sum([np.cos(2*(theta_i-alpha))**2 for theta_i in Bi], 
                   axis = 0)
    
    pB = Numerator/Denom
    
    return pB

#Equation 7
def pB_MZP_sum(Bi, alpha):
    '''Converts Solar position angle of an image point,`alpha`, into Coronal 
    polarized brightness, `pB`. Uses summation as another way to calculate cornal polarized brightness.
    
    This function takes in the var of `alpha`.
    
    Parameters
    ----------
    alpha : np.ndarray
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.
    Raises
    none
    ------
    
    '''
    
    
    pB = (-4/3)*(np.sum([ith_polarizer_brightness 
                         * np.cos(2*(ith_angle - alpha)) 
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    
    return pB


#Equation 8
def polarized_B(pB, B_theta, theta, alpha):
    '''Converts Radiance through a polarizer at angle theta,`B_theta`,
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
    Raises
    none
    ------
    
    '''
    
    B = 2*B_theta + (np.cos(2*(theta - alpha)))*pB
    
    return B


#Equation 9
def B_MZP_sum(Bi):
    '''Converts unpolarized brightness,`B`, into an average unpolarized brightness. 
    Uses summation averaging over the three polarizer positions M,Z, and P.
    
    This function takes in a var of, `Bi`, which is a dictionary of these different polarizer positions.
    
    Parameters
    ----------
    Bi : dict
    
    Returns
    -------
     float
        The float that is returned is defined to be var, `pB`.
    Raises
    none 
    ------
    
    '''
    
    B = (2/3)*(np.sum([ith_polarizer_brightness  
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    
    return B

#Equation 10
def pB_prime(Bi, alpha):
    '''Converts total polarizer brightness at an ith angle, `Bi`, and solar position angle of an image point,
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
    Raises
    ? 
    ------
    
    '''
    
    pB_prime = (-4/3)*(np.sum([ith_polarizer_brightness 
                         * np.sin(2*(ith_angle - alpha)) 
                        for ith_angle, 
                         ith_polarizer_brightness in Bi.items()], axis=0))
    return pB_prime


#Equation 11
def B_theta_using_pB_prime(B, pB, pB_prime, theta, alpha):
    '''Converts unpolarized brightness,`B`, original polarized brightnesses,`pB`, `pB_prime`, 
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
        
    Raises
    ValueError if shapes do match
    ------
    
    '''
    #checks to see if all np.ndarrarys match each other in shape
    # if B.shape != theta.shape:
    #     raise ValueError(f'B and theta shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'theta.shape = {theta.shape}')
        
    # elif B.shape != alpha.shape:
    #     raise ValueError(f'B and alpha shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif B.shape != pB.shape:
    #     raise ValueError(f'B and pB shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'pB.shape = {pB.shape}')
        
    # elif B.shape != pB_aster.shape:
    #     raise ValueError(f'B and pB_aster shapes do not match. ' 
    #                      f'B.shape = {B.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
        
    # elif theta.shape != alpha.shape:
    #     raise ValueError(f'theta and alpha shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif theta.shape != pB.shape:
    #     raise ValueError(f'theta and pB shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'pB.shape = {pB.shape}')
        
    # elif theta.shape != pB_aster.shape:
    #     raise ValueError(f'theta and pB_aster shapes do not match. ' 
    #                      f'theta.shape = {theta.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
        
    # elif alpha.shape != pB.shape:
    #     raise ValueError(f'alpha and pB shapes do not match. ' 
    #                      f'alpha.shape = {alpha.shape} '
    #                      f'pB_aster.shape = {pB.shape}')
        
    # elif alpha.shape != pB_aster.shape:
    #     raise ValueError(f'alpha and pB_aster shapes do not match. ' 
    #                      f'alpha.shape = {alpha.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
    # elif pB.shape != pB_aster.shape:
    #     raise ValueError(f'pB and pB_aster shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
    
    B_theta = (1/2)*(B - np.cos(2*(theta - alpha))*pB 
                    - np.sin(2*(theta - alpha))*pB_prime)
    
    return B_theta


#Equation 12
def Stokes_Q(B_z,B_m,B_p):
    '''Converts unpolarized brightness zero,`B_z`, unpolarized brightness plus direction,`B_p`,
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
        
    Raises
    ValueError if shapes do match
    ------
    
    '''
    
    #checks to see if all np.ndarrarys match each other in shape
    # if B_z.shape != B_m.shape:
    #     raise ValueError(f'B_z and B_m shapes do not match. ' 
    #                      f'B_z.shape = {B_z.shape} '
    #                      f'B_m.shape = {B_m.shape}')
        
    # elif B_z.shape != B_p.shape:
    #     raise ValueError(f'B_z and B_p shapes do not match. ' 
    #                      f'B_z.shape = {B_z.shape} '
    #                      f'B_p.shape = {B_p.shape}')
        
    # elif B_m.shape != B_p.shape:
    #     raise ValueError(f'B_m and B_p shapes do not match. ' 
    #                      f'B_m.shape = {B_m.shape} '
    #                      f'B_p.shape = {B_p.shape}')
        
    Q = (2/3)*(2*B_z - B_m - B_p)
    
    return Q


#Equation 13
def Stokes_U(B_p,B_m):
    '''Converts unpolarized brightness in the plus direction,`B_p`, unpolarized brightness in the minus direction, `B_p`
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
        
    Raises
    ValueError if shapes do match
    ------
    
    '''
    
    #checks to see if all np.ndarrarys match each other in shape
    # if B_m.shape != B_p.shape:
    #     raise ValueError(f'B_m and B_p shapes do not match. ' 
    #                      f'B_m.shape = {B_m.shape} '
    #                      f'B_p.shape = {B_p.shape}')
    
    U = (2/np.sqrt(3))*(B_p - B_m)
    
    return U


#Equation 15
def theta_maximum(pB, pB_prime, alpha):
    '''Converts Coronal polarized brightness, `pB`, `pb_prime`, 
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
        
    Raises
    ValueError if shapes do match
    ------
    
    '''
    
    #checks to see if all np.ndarrarys match each other in shape
    # if pB.shape != pB_aster.shape:
    #     raise ValueError(f'pB and pB_aster shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
        
    # elif pB.shape != alpha.shape:
    #     raise ValueError(f'pB and alpha shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    # elif pB_aster.shape != alpha.shape:
    #     raise ValueError(f'pB_aster and alpha shapes do not match. ' 
    #                      f'pB_aster.shape = {pB_aster.shape} '
    #                      f'alpha.shape = {alpha.shape}')
        
    theta_max = ((1/2)*np.arctan(pB_prime/pB)) + np.pi/2 + alpha
    
    return theta_max


#Equation 16
def polarization_fraction(B, pB, pB_prime):
    '''Converts Coronal polarized brightness, `pB`,`pB_prime`, and
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
        
    Raises
    ValueError if shapes do match
    ------
    
    '''
    
    #checks to see if all np.ndarrarys match each other in shape
    # if pB.shape != pB_aster.shape:
    #     raise ValueError(f'pB and pB_aster shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'pB_aster.shape = {pB_aster.shape}')
        
    # elif pB.shape != B.shape:
    #     raise ValueError(f'pB and B shapes do not match. ' 
    #                      f'pB.shape = {pB.shape} '
    #                      f'B.shape = {B.shape}')
        
    # elif pB_aster.shape != B.shape:
    #     raise ValueError(f'pB_aster and B shapes do not match. ' 
    #                      f'pB_aster.shape = {pB_aster.shape} '
    #                      f'B.shape = {B.shape}')
        
    p = np.sqrt((pB**2 + pB_prime**2)) / B
    
    return p

 
#TODO think about odd array lengths (rounding error).
def alpha_func(array_length):
    '''Makes an array of alpha values on a mesh grid. This function take in the parameter of 
    arrary_length to properly parameterize the created meshgrid.
    
    This function takes in a var of `array_length`.
    
    Parameters
    ----------
    array_length : int 
    
    
    Returns
    -------
     np.array
        The np.array that is returned is defined to be var, `alpha`.
        
    Raises
    nothing
    ------
    
    '''
    x = np.arange(-array_length//2,array_length//2)
    y = np.arange(-array_length//2,array_length//2)
    xx, yy = np.meshgrid(x, y)
    alpha = np.rot90(np.arctan2(yy,xx), k = 1)
    alpha=np.fliplr(alpha)
    return alpha
