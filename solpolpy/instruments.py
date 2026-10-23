"""Instrument specific code."""

from os import PathLike

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.wcs import WCS
from ndcube import NDCollection, NDCube

from solpolpy.errors import TooFewFilesError, UnsupportedInstrumentError


def load_punch(input_data: str | PathLike | NDCube | list, mask: np.ndarray | None = None) -> NDCollection:
    """Load one PUNCH mosaic or a supplied M/Z/P triplet into solpolpy.

    Parameters
    ----------
    input_data : path-like, ndcube.NDCube, or sequence
        A PTM, PIM, or PAM FITS path or already-loaded cube, or three
        individual M/Z/P paths or cubes in any order. Punchbowl
        ``PUNCHCube`` objects and their metadata fields are also accepted.
        Supported standard layouts:
        1. Level 1 separate M/Z/P frames,
        2. Level 2 PTM with M/Z/P layers,
        3. Level 3 PIM with M/Z/P layers, and
        4. Level 3 PTM/PAM with tB/pB/pBp layers (pBp denotes pB-prime).
    mask : numpy.ndarray, optional
        Additional 2D boolean mask shared by all channels;
        True marks invalid pixels. Combined with existing and invalid-data masks.

    Returns
    -------
    ndcube.NDCollection
        Three 2D cubes keyed by M/Z/P or tB/pB/pBp, ready for ``solpolpy.resolve``.
    """
    if isinstance(input_data, (str, PathLike, NDCube)):
        input_data = [input_data]
    else:
        input_data = list(input_data)
    if len(input_data) not in (1, 3):
        raise ValueError("Supply one PUNCH mosaic or three individual M/Z/P frames")

    channels = {}
    angles = {"M": -60, "Z": 0, "P": 60}
    for source in input_data:
        if isinstance(source, (str, PathLike)):
            with fits.open(source) as hdul:
                header = hdul[1].header.copy()
                data = hdul[1].data.copy()
                # Passing the HDU list retains any spatial distortion tables.
                wcs = WCS(header, hdul, naxis=2)
                sigma = None
                if "Uncertainty array" in hdul:
                    ratio = hdul["Uncertainty array"].data
                    if ratio.shape != data.shape:
                        raise ValueError("PUNCH uncertainty array must match the signal shape")
                    sigma = np.full(data.shape, np.inf, dtype=float)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        np.divide(data, ratio, out=sigma, where=ratio != 0)
                source_mask = None
        elif isinstance(source, NDCube):
            # NormalizedMetadata exposes MetaField objects rather than raw values.
            header = {}
            for key in source.meta:
                value = source.meta[key]
                header[key] = value if isinstance(value, u.Quantity) else getattr(value, "value", value)
            data, wcs = source.data, source.wcs
            if wcs.pixel_n_dim == 3:
                wcs = wcs.dropaxis(2)
            sigma = None
            if source.uncertainty is not None:
                sigma = (
                    source.uncertainty.array
                    if isinstance(source.uncertainty, StdDevUncertainty)
                    else source.uncertainty.represent_as(StdDevUncertainty).array
                )
            source_mask = source.mask
        else:
            raise TypeError("PUNCH inputs must be FITS paths or NDCube objects")

        if data.ndim == 3 and data.shape[0] == 3 and len(input_data) == 1:
            labels = [str(header.get(f"OBSLAYR{i}", "")).removeprefix("Polar_") for i in (1, 2, 3)]
            if set(labels) not in ({"M", "Z", "P"}, {"tB", "pB", "pBp"}):
                raise ValueError("Mosaic OBSLAYR1-3 must identify M/Z/P or tB/pB/pBp planes")
        elif data.ndim == 2 and len(input_data) == 3:
            angle = u.Quantity(header["POLAR"], u.deg).to_value(u.deg)
            labels = [key for key, value in angles.items() if np.isclose(angle, value, rtol=0, atol=1e-6)]
            if len(labels) != 1:
                raise ValueError("Individual PUNCH frames must have POLAR = -60, 0, or 60 degrees")
        else:
            raise ValueError("Expected one (3, ny, nx) mosaic or three (ny, nx) frames")

        if wcs.pixel_n_dim != 2:
            raise ValueError("PUNCH channels require a two-dimensional spatial WCS")
        if mask is not None and (np.shape(mask) != data.shape[-2:] or np.asarray(mask).dtype != bool):
            raise ValueError("mask must be a boolean array matching the spatial image shape")

        for i, label in enumerate(labels):
            if label in channels:
                raise ValueError(f"Duplicate PUNCH channel: {label}")
            signal = data[i] if data.ndim == 3 else data
            error = None if sigma is None else (sigma[i] if data.ndim == 3 else sigma)
            invalid = ~np.isfinite(signal)
            if source_mask is not None:
                invalid |= np.broadcast_to(source_mask, data.shape)[i] if data.ndim == 3 else source_mask
            if error is not None:
                invalid |= ~np.isfinite(error) | (error < 0)
            if mask is not None:
                invalid |= mask
            meta = dict(header)
            if label in angles:
                if str(meta.get("POLARREF", "")).lower() not in ("solar", "instrument"):
                    raise ValueError("M/Z/P inputs require POLARREF = Solar or Instrument")
                meta["POLAR"] = angles[label] * u.deg
            else:
                meta["POLAR"] = label
            channels[label] = NDCube(
                signal.copy(),
                wcs=wcs.deepcopy(),
                meta=meta,
                mask=invalid,
                uncertainty=None if error is None else StdDevUncertainty(error.copy()),
            )

    order = ("M", "Z", "P") if "M" in channels else ("tB", "pB", "pBp")
    if "M" in channels and len({str(cube.meta["POLARREF"]).lower() for cube in channels.values()}) != 1:
        raise ValueError("M/Z/P frames must share the same POLARREF")
    return NDCollection([(key, channels[key]) for key in order], meta={}, aligned_axes="all")


def get_data_angle(header):
    angle_hdr = header["POLAR"]
    if isinstance(angle_hdr, float):
        angle = angle_hdr * u.degree
    elif isinstance(angle_hdr, str):
        angle = float(angle_hdr.split("Deg")[0]) * u.degree
    else:
        try:
            angle = float(angle_hdr) * u.degree
        except ValueError:
            msg = "Polar angle in the header could not be read for this instrument."
            raise UnsupportedInstrumentError(msg)

    return angle


def load_data(
    path_list: list[str], mask: np.ndarray | None = None, use_instrument_mask: bool = False, hdu_index: int = 0
) -> NDCollection:
    """Basic loading function. See `load_with_occulter_mask`.

    Parameters
    ----------
    path_list: List[str]
        list of paths to FITS files to be loaded

    mask: Optional[np.ndarray]
        An optional mask to couple with images. Overrides any mask computed when use_instrument_mask=True.

    use_instrument_mask: bool
        If true, loads an instrument mask for common instruments defined in `get_instrument_mask`.

    hdu_index: int
        The index of the HDU (zero-based) to load data from.
        For many spacecraft this should be 0. For PUNCH it should be set to 1.

    Returns
    -------
    NDCollection
        The data are loaded as NDCollection object with WCS and header information available.
        The keys are labeled as 'angle_1', 'angle_2, 'angle_3', ...

    """
    # get length of list to determine how many files to process.
    if len(path_list) < 2:
        msg = "Requires at least 2 FITS files"
        raise TooFewFilesError(msg)

    data_out = []
    for _i, data_path in enumerate(path_list):
        with fits.open(data_path) as hdul:
            wcs = WCS(hdul[hdu_index].header)

            if use_instrument_mask and mask is None:
                mask = get_instrument_mask(hdul[hdu_index].header)

            if mask is None:  # make a mask of False if none is provided
                mask = np.zeros(hdul[hdu_index].data.shape, dtype=bool)

            angle = get_data_angle(hdul[hdu_index].header)

            data_out.append(
                (str(angle), NDCube(hdul[hdu_index].data, mask=mask, wcs=wcs, meta=dict(hdul[hdu_index].header)))
            )

    return NDCollection(data_out, meta={}, aligned_axes="all")


def construct_mask(
    inner_radius: float, outer_radius: float, center_x: int, center_y: int, shape: tuple[int, int]
) -> np.ndarray:
    """Constructs a mask where False indicates the pixel is valid and True masks the pixel as invalid.

    Pixels with radial distance between `inner_radius` and `outer_radius` are valid. Every other pixel is invalid.
    This is a standard coronograph mask.

    Parameters
    ----------
    inner_radius : float
        inner radius of the mask. pixels closer to the center than this are masked as invalid
    outer_radius : float
        outer radius of the mask. pixels farther from the center than this are masked as invalid
    center_x : int
        center pixel x coordinate
    center_y : int
        center pixel y coordinate
    shape : Tuple[int, int]
        the image dimensions

    Returns
    -------
    np.ndarray
        a coronograph mask that marks invalid pixels
        (those with radius less than `inner_radius` or greater than `outer_radius`) as True

    """
    xx, yy = np.ogrid[0 : shape[0], 0 : shape[1]]
    mask = np.zeros(shape, dtype=bool)
    mask[(xx - center_x) ** 2 + (yy - center_y) ** 2 < inner_radius] = True
    mask[(xx - center_x) ** 2 + (yy - center_y) ** 2 > outer_radius] = True
    return mask


def get_instrument_mask(header: fits.Header) -> np.ndarray:
    """Gets a coronograph mask for common instruments.

    Supports common solar instruments: LASCO, COSMO K, SECCHI

    Parameters
    ----------
    header : astropy.io.fits.Header
        a FITS header for the file in question

    Returns
    -------
    np.ndarray
        a pixel mask where valid pixels are marked False, masked out invalid pixels are flagged as True

    """
    x_shape = header["NAXIS1"]
    y_shape = header["NAXIS2"]

    center_x = header["CRPIX1"]
    center_y = header["CRPIX2"]

    R_sun = 960.0 / header["CDELT1"]  # TODO: clarify where 960. comes from

    radius = 0
    if header["INSTRUME"] == "LASCO":
        if "C2" in header["DETECTOR"]:
            radius = 2.5
        elif "C3" in header["DETECTOR"]:
            radius = 4.0
        else:
            msg = "The data in this file is not formatted according to a known instrument."
            raise UnsupportedInstrumentError(msg)
    elif header["INSTRUME"] == "COSMO K-Coronagraph":
        radius = 1.15
    elif header["INSTRUME"] == "SECCHI":
        if "COR1" in header["DETECTOR"]:
            radius = 1.57
        elif "COR2" in header["DETECTOR"]:
            radius = 3.0
        else:
            msg = "The data in this file is not formatted according to a known instrument."
            raise UnsupportedInstrumentError(msg)
    else:
        msg = "The data in this file is not formatted according to a known instrument."
        raise UnsupportedInstrumentError(msg)

    y_array = [(j_step - center_y) / R_sun for j_step in range(y_shape)]
    outer_edge = np.max(y_array) * R_sun

    return construct_mask((radius * R_sun) ** 2, outer_edge**2, center_x, center_y, (x_shape, y_shape))
