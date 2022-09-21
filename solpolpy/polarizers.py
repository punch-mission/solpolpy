#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solpolpy subclasses for basic coordinate transforms
"""

import numpy as np


class Br(Solpolpy):
    '''
    solpolpy.Br - Converts unpolarized brightness,`B`, and Coronal polarized brightness,
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
    ValueError
      Dimensional mismatch or non-array inputs cause this to be thrown.
          
    Notes
    ------
    Equation 1 in Deforest et al. 2022.
    """
    '''
    def __init__(self, data: np.ndarray):
        (B - pB) / 2
    
    def apply( self, data: np.ndarray ):
        return( data )

    def __str__(self):
        self._strtmp = f"Br ({self.params['Br']})"
        return super().__str__()