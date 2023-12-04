import warnings
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from ndcube import NDCube
from astropy.wcs import WCS
from astropy.io import fits
from sunpy.map.sources import LASCOMap



def create_oculter_mask(file_path, 
                        mask_radii=None, 
                        color_map=None, 
                        centerpx_x=None,
                        centerpx_y=None,
                        solar_radii=None):
    '''
    create_oculter_mask creates a simple masked dataset for coronagraph data, 
    it applies masks and color maps.

    Input
    -----

    file_path : path to a FITS file to be masked

    mask_radii : specifies the size of the mask in Solar Radii. If input overwrites
        default mask sizes dependent on instruments.

    color_map : specifies color maps to be applied to data

    centerpx_x : if defined specifies the x position center of the masked region

    centerpx_y : if defined specifies the y position center of the masked region

    solar_radii : specify a solar radii in pixels


    Output
    ------

    output_ndcube : returns an ndcube of output data, wcs and mask

    '''

    file_of_interest = fits.open(file_path)
    file_data = file_of_interest[0].data
    file_meta = file_of_interest[0].header

    if file_meta['CRPIX1']:
        center_x = file_meta['CRPIX1']
    else:
        center_x=int(file_data.shape[0]/2)
        warnings.warn("No CRPIX1 found in FITS header, setting CRPIX1 to center pixel")
    
    if file_meta['CRPIX2']:
        center_y = file_meta['CRPIX2']
    else:
        center_y=int(file_data.shape[1]/2)
        warnings.warn("No CRPIX2 found in FITS header, setting CRPIX2 to center pixel")

    if file_meta['CDELT1']:
        R_sun = 960. / file_meta['CDELT1'] # this needs to be fixed with explanation, always the same?
    else: 
        R_sun=50
        warnings.warn("No CDELT1 found in FITS header, setting solar_radii to 50")

    if file_meta['INSTRUME'] == 'LASCO':
        detector_name = file_meta['DETECTOR']
        if 'C2' in detector_name:
            mask_R = 2.5
            
        elif 'C3' in detector_name:
            mask_R = 4.0

    elif file_meta['INSTRUME'] == 'COSMO K-Coronagraph':
        mask_R = 1.15

    elif file_meta['INSTRUME'] == 'SECCHI':
        detector_name = file_meta['DETECTOR']
        if 'COR1' in detector_name:
            mask_R = 1.57

        elif 'COR2' in detector_name:
            mask_R = 3.0

    else:
        warnings.warn("No INSTRUME found in FITS header, setting mask_radii to 1")
        mask_R=1.0


    # user defined inputs to overwrite pre-defined inputs
    if mask_radii != None:
        mask_R=mask_radii

    if centerpx_x != None:
        center_x=centerpx_x

    if centerpx_y != None:
        center_y=centerpx_y


    # obtain array dimensions
    x_shape, y_shape = file_data.shape

    # Calibrate to solar radius
    x_array = np.empty(x_shape)
    y_array = np.empty(y_shape)

    for i_step in range(x_shape):
        x_array[i_step] = (i_step - center_x) / R_sun
        
    for j_step in range(y_shape):
        y_array[j_step] = (j_step - center_y) / R_sun

    
    rr = np.max(y_array) * R_sun
    xx, yy = np.ogrid[0:x_shape, 0:y_shape]
    output_mask = np.ones(file_data.shape)  # the final mask array
    mask = (xx - center_x) ** 2 + (yy - center_y) ** 2 < (mask_R * R_sun) ** 2
    output_mask[mask] = 0
    mask = (xx - center_x) ** 2 + (yy - center_y) ** 2 > rr ** 2
    output_mask[mask] = 0


    output_ndcube = NDCube(file_data, WCS(file_meta))
    output_ndcube.mask = output_mask
    output_ndcube.meta = file_meta


    return output_ndcube




def create_color_map_name(file_path):
    '''
    A simple wrapper to create a color map name from an input FITS file

    Input
    -----

    file_path : path to a FITS file to be masked


    Output
    ------

    color_map : returns the colormap string

    '''
    
    file_of_interest = fits.open(file_path)
    file_data = file_of_interest[0].data
    file_meta = file_of_interest[0].header


    if file_meta['INSTRUME'] == 'LASCO':
        detector_name = file_meta['DETECTOR']
        if 'C2' in detector_name:
            color_map='soholasco2'
            
        elif 'C3' in detector_name:
            color_map='soholasco3'

    elif file_meta['INSTRUME'] == 'COSMO K-Coronagraph':
        color_map='kcor'

    elif file_meta['INSTRUME'] == 'SECCHI':
        detector_name = file_meta['DETECTOR']
        if 'COR1' in detector_name:
            color_map='stereocor1'

        elif 'COR2' in detector_name:
            color_map='stereocor2'

    else:
        warnings.warn("No INSTRUME found in FITS header, setting color_map soholasco2")
        color_map='soholasco2'

    output_color_map = plt.get_cmap(color_map)

    return output_color_map
