import warnings

from astropy.io import fits


def get_colormap_str(file_path: str) -> str:
    """Retrieve a color map name from an input FITS file
    Parameters
    ----------
    file_path : str
        path of file to open

    Returns
    -------
    str
        name of appropriate colormap
    """
    file_of_interest = fits.open(file_path)
    file_meta = file_of_interest[0].header

    if file_meta['INSTRUME'] == 'LASCO':
        detector_name = file_meta['DETECTOR']
        if 'C2' in detector_name:
            color_map = 'soholasco2'
        elif 'C3' in detector_name:
            color_map = 'soholasco3'
        else:
            warnings.warn("No valid instrument found, setting color_map soholasco2")
            color_map = 'soholasco2'
    elif file_meta['INSTRUME'] == 'COSMO K-Coronagraph':
        color_map = 'kcor'
    elif file_meta['INSTRUME'] == 'SECCHI':
        detector_name = file_meta['DETECTOR']
        if 'COR1' in detector_name:
            color_map = 'stereocor1'
        elif 'COR2' in detector_name:
            color_map = 'stereocor2'
        else:
            warnings.warn("No valid instrument found, setting color_map soholasco2")
            color_map = 'soholasco2'
    else:
        warnings.warn("No valid instrument found, setting color_map soholasco2")
        color_map = 'soholasco2'

    return color_map
