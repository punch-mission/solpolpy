from astropy.io import fits
import numpy


def resolve(input_data, out_polarize_state, separation=None, alpha=None, Error=False):
    """
    Apply - apply a depolarization transform to a set of input
    dataframes.

    Parameters
    ----------
    input_data : Dict[str, np.ndarray]
        Dictionary formated as follows:

        - Stokes dictionary
        "I":np.array - Should be included as a triplet of I,Q,U, and optionally V
        "Q":np.array - Should be included as a triplet of I,Q,U, and optionally V
        "U":np.array - Should be included as a triplet of I,Q,U, and optionally V
        "V":np.array - Should be included as a triplet of I,Q,U, and optionally V

        - Brightness & Polarized brightness dictionary
        "B":np.array - Should be included as a double of B, pB
        "pB":np.array - Should be included as a doublet of B, pB

        - Radial & tangential brightness dictionary
        "Br":np.array - Should be included as a doublet of Bt, Br
        "Bt":np.array - Should be included as a doublet of Bt, Br

        - MZP triplet dictionary
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

    out_polarize_state : string
      This is the polarization state you want to convert your input
      dataframes to.  These include:

    - Stokes: convert the input dataframes to output Stokes dataframes. If
        3 polarization angles are input, or can be derived, the I, Q, and U
        Stokes parameters are output. If 4 or more polarization angles are
        provided, and the angles are conducsive, the full I,Q,U, and V
        Stokes parameters are provided.

    - B: converts input dataframes into B ("unpolarized
        brightness") parameters.

    - pB: converts input dataframes into pB ("polarized brightness")
        parameters.

    - Bt: Produces "Bt", the radiance observed through a linear polarizer
        oriented tangentially to a solar-concentric circle passing through
        the image point of interest.

    - Br: Produces "Br" the radiance observed through a linear polarizer
        oriented radially to the centre Sun through the same point.

    - BrBt: Produces both the  radiance observed through a linear polarizer
        oriented tangentially to a solar-concentric circle passing through
        the image point of interest "Bt", and the radiance observed through
        a linear polarizer oriented radially to the centre Sun through the
        same point "Br".

    - 3pol or MZP: converts input dataframes into a system of virtual
        polarizer triplets each separated by 60 degrees (Minus [-60],
        Zero [0], Plus [60]).

    - 4pol: converts input dataframes into a system of virtual polarizer
        triplets each separated by 45 degrees at -45, 0, 45, 90.

    - Xpol: converts input dataframes into a system of virtual polarizer
        triplets each separated by X degrees, specified by the optional
        input separation, which should be input in degrees.

    - BpB: converts input dataframes into two dataframes, B ("unpolarized
        brightness") parameter and its counterpart pB ("polarized brightness")

    Raises
    ------
    AssertionError
      This gets raised if the data cannot be converted or polarization
      transformation cannot calculated due to a discontinuity or infinity.

    Returns
    -------
    numpy.ndarray
        The transformed vector data are returned as a numpy.ndarray.  Most
        Transforms maintain the dimensionality of the source vectors.  Some
        embed (increase dimensionality of the vectors) or project (decrease
        dimensionality of the vectors); additional input dimensions, if
        present, are still appended to the output vectors in all any case.
    """
    pass

    # # create output dictionary
    # #data_out={}
    #
    #
    # # check if go a dictionary
    # print(type(data_in))
    # if isinstance(data_in, dict):
    #     data_out=data_in
    #     print('got dict')
    #
    #
    # # check if got a list
    # if isinstance(data_in, list):
    #     data_out={}
    #     list_len=len(data_in)
    #     assert list_len >= 2, 'requires at least 2 FITS files'
    #
    #     for xlist_item in data_in:
    #         with fits.open(xlist_item) as hdul:
    #             assert hdul[0].header['INSTRUME'] == 'SECCHI', 'requires FITS to be SECCHI COR data files'
    #             data_out[hdul[0].header['POLAR']]=hdul[0].data
    #             image_hdr = hdul[0].header
    #     print('got list')
    #
    #
    # return data_out
