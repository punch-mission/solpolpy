from typing import List

from ndcube import NDCube, NDCollection
from astropy.io import fits
from astropy.wcs import WCS


def load_data(path_list: List[str]) -> NDCollection:
    """
    Parameters
    ----------
    path_list: List[str] 
        list of paths to be loaded

    Returns
    -------
    NDCollection
        The data are loaded as NDCollection object with WCS and header information available.
        The keys are labeled as 'angle_1', 'angle_2, 'angle_3', ...
    """
    # create list of FITS
    fits_type = []

    # get length of list to determine how many files to process.
    list_len = len(path_list)
    assert list_len >= 2, 'requires at least 2 FITS files'
    
    data_out = []
    for i, data_path in enumerate(path_list):
        with fits.open(data_path) as hdul:
            wcs = WCS(hdul[0].header)
            data_out.append(("angle_" + str(i), 
                             NDCube(hdul[0].data, 
                                    wcs=wcs, 
                                    meta=hdul[0].header)))
            
    return NDCollection(data_out, meta={}, aligned_axes="all")
