#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solpolpy subclasses for basic coordinate transforms
"""

from astropy.io import fits
import numpy

#TODO: so this should take:

#Instead of dict to feed in, we create data object whiych gets passed into the depolariser object

#Counts number of arrays. And decides if IQUV or IQU
#Stokes  data cube
#A couple of common basis
#Separation could be list or numeric value for angles
#Alpha=keyword or array or number

#Float or int create as a numpy array or FITS or WCS should be passed in as dict

#Core - core functionality - calls Polarizers
#Polarizers - collection of conversions
#Errors - collection of Error conversions
#Validation - checks the dictionary is fine and what can be generated.

#Apply 



def data_type(*datain):
    '''data_type - converts several input datatypes into a dictionary format 
    passable into the Depolarize object.
    
    The Depolarize object requires specific information with data frames 
    to be depolarized, such as polarizer angles, or data types.  The provided
    data will be passed as a dictionary with angle or data type key and data 
    pairs. However data may be available in many formats, such as FITS.
    
    The output is a key:dataframe dictionary. dataframes, will be in numpy 
    array format.  The different key pairs will be:

    - dictionary
    "I":np.array - Should be included as a triplet of I,Q,U, and optionally V
    "Q":np.array - Should be included as a triplet of I,Q,U, and optionally V
    "U":np.array - Should be included as a triplet of I,Q,U, and optionally V
    "V":np.array - Should be included as a triplet of I,Q,U, and optionally V

    - Stokes dictionary
    "Stokes":np.array[x3] - the array will contain three data frames, I, Q, U in that order
    "Stokes":np.array[x4] - the array will contain four data frames, I, Q, U, and V in that order

    - Brightness & Polarized brightness dictionary
    "B":np.array - Should be included as a double of B, pB
    "pB":np.array - Should be included as a doublet of B, pB

    - Combined brightness & polarized brightness dictionary
    "BpB":np.array[x2] - the array will contain two data frames, B and pB,
        in that order

    - Radial & tangential brightness dictionary
    "Br":np.array - Should be included as a doublet of Bt, Br
    "Bt":np.array - Should be included as a doublet of Bt, Br

    - Combined radial & tangential brightness dictionary
    "BrBt":np.array[x2] - the array will contain two data frames, Br and Bt
        in that order

    - MSP triplet dictionary
    "M":np.array - Should be included as a triplet of M,Z,P
    "Z":np.array - Should be included as a triplet of M,Z,P
    "P":np.array - Should be included as a triplet of M,Z,P

    - Angular dictionary [where keys are ints or floats representing 
        polarizer angle]
    X1:np.array[xn] - Should be included as at least a triplet of angles
    X2:np.array[xn] - Should be included as at least a triplet of angles
    X3:np.array[xn] - Should be included as at least a triplet of angles
         ...                               ...
    Xn:np.array[xn] - Should be included as at least a triplet of angles

    - Combined MZP triplet dictionary
    "MZP":np.array[x3] - the array will contain three data frames, M,Z,and P,
        in that order

    - Three Polarizer dictionary
    "3Pol":np.array[x3] - the array will contain three data frames separated 
                          by 60 degress [-60,0,60], in that order
    - dictionary
    "4Pol":np.array[x4] - the array will contain frames data frames separated 
                          by 45 degress [-45,0,45,90], in that order

     
    Examples
    --------


    
    Methods
    -------
        
         
    Built-in Subclasses
    -------------------
    
    The Depolarize class is a container for the subclasses that do the actual 
    mathematical work.  Each subclass represents a depolarization operation; 
    an instance of that subclass represents a single depolarization selected 
    from the family by the depolarizers you pass to the constructor.
    
    The depolarizer module itself defines the following subclasses:
                    
        - Composition: groups together multiple depolarization transforms 
        into one composite depolarization transform.
          
        
    '''

    count=0
    for sep_datain in datain:

        if count==0:
            if isinstance(sep_datain, fits.hdu.hdulist.HDUList):
                print("got fits", count)
        if count>0:
            if isinstance(sep_datain, fits.hdu.hdulist.HDUList):
                print("got fits", count)
        count=count+1

        pass
    


    dataout=datain

    return datain

