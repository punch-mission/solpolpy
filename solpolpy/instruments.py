
from typing import List
from ndcube import NDCube, NDCollection
from astropy.io import fits
from astropy.wcs import WCS


def load_data(path_list: List[str]) -> NDCollection:
    """
    path_list: String of path list where group of files to be loaded is present.
    default_alpha: bool, optional
                   If default_alpha is true the alpha matrix will be referenced from Solar North. Default is True.
                   Make it False to provide alpha (in FITS format) with same size as of input data
                   or select the alpha matrix (from the list).

    Returns
    -------
    NDCollection
        The data are loaded as NDCollection object with WCS and header information available.
        The keys are labelled as 'angle_1', 'angle_2, 'angle_3', ...
        Alpha matrix is also added apart from the polarizing angles data and can be accessed through 'alpha' key.
    """

    # create list of FITS
    fits_type = []

    # get length of list to determine how many files to process.
    list_len = len(path_list)
    assert list_len >= 2, 'requires at least 2 FITS files'

    # for xlist_item in path_list:
    #     with fits.open(xlist_item) as hdul:
    #         fits_type.append(hdul[0].header['DETECTOR'])
    #
    # if len(set(fits_type)) != 1:
    #     raise Exception("Input FITS are of different types")

    data_out = []
    i = 0

    for data_path in path_list:
        with fits.open(data_path) as hdul:
            i = i+1
            wcs = WCS(hdul[0].header)
            data_out.append(("angle_" + str(i), NDCube(hdul[0].data, wcs=wcs, meta=hdul[0].header)))
    #         size = hdul[0].data.shape
    # if default_alpha:
    #     alph = radial_north(size)
    #     data_out.append(("alpha", NDCube(alph, wcs=wcs)))
    # else:
    #     inp = input("Do you wish to provide an alpha array?").lower()
    #     if inp.startswith('n'):
    #         print("Continuing with default options... Waiting for input...")
    #         inp_ref = input("Choose the reference along the Solar: \"North\" or \"West\":").lower()
    #         if inp_ref.startswith('n'):
    #             alpha = radial_north(size)
    #         else:
    #             alpha = radial_west(size)
    #         data_out.append(("alpha", NDCube(alpha, wcs=wcs, meta=hdul[0].header)))
    #     elif inp.startswith('y'):
    #         print("Provide the alpha matrix in FITS format")
    #         alpha_path = input("Provide the path of alpha FITS file:")
    #         hdu = fits.open(alpha_path)
    #         if np.max(hdu[0].data) > 2 * np.pi + 1:
    #             alph = hdu[0].data * deg2rad
    #         else:
    #             alph = hdu[0].data
    #         data_out.append(("alpha", NDCube(alph, wcs=wcs, meta=hdu[0].header)))

    print("Hurray!!! Data loaded successfully.")
    return NDCollection(data_out, meta={}, aligned_axes="all")
