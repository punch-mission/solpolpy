import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import sunpy.visualization.colormaps as cm  # noqa: F401
from astropy.io import fits
from ndcube import NDCollection, NDCube
from sunkit_image.radial import fnrgf, intensity_enhance, nrgf, rhef
from sunkit_image.utils import equally_spaced_bins


def plot_collection(collection,
                    figsize=(8, 8),
                    show_colorbar=False,
                    lat_ticks=None,
                    lon_ticks=None,
                    major_formatter="dd",
                    xlabel="HP Longitude",
                    ylabel="HP Latitude",
                    vmin=None,
                    vmax=None,
                    cmap="Greys_r",
                    ignore_alpha=True,
                    fontsize=18,
                    **kwargs):
    """Plot a solpolpy NDCollection input or output.

    Parameters
    ----------
    collection : NDCollection, ndarray, or 3D color_image.
        collection to visualize
    figsize : Tuple[float, float]
        figure size according to Matplotlib
    show_colorbar : bool
        whether to show a colorbar
    lat_ticks : Optional[np.ndarray]
        if provided, shows as the tick marks for latitude. default values used otherwise.
    lon_ticks : Optional[np.ndarray]
        if provided, shows as the tick marks for longitude. default values used otherwise.
    major_formatter : str
        the formatter for major tickmarks as specified by Matplotlib
    xlabel : str
        label for plot x axes
    ylabel : str
        label for plot y axes
    vmin : float, list of floats, or None
        minimum values of the plots. if a list is provided, they are applied left to right to each plot
    vmax : float, list of floats, or None
        maximum values of the plots. if a list is provided, they are applied left to right to each plot
    cmap : str or Matplotlib colormap
        a Matplotlib accepted colormap or colormap string
    ignore_alpha : bool
        whether to plot the alpha array. defaults to True as it is not normally helpful to visualize.
    fontsize : int
        font size for some aspects of the plot
    **kwargs : Additional imshow keyword arguments.
            Extra parameters to pass to `imshow()`.

    Returns
    -------
    Matplotlib figure and axes
        the plotted figure and axes are returned for any additional edits

    """
    # Check if collection is an NDCollection, ndarray, or color_image
    if isinstance(collection, dict):
        collection_keys = list(collection.keys())
        if ignore_alpha:
            collection_keys = [k for k in collection_keys if k != "alpha"]
        ax_count = len(collection_keys)
        first_item = collection[collection_keys[0]]
        wcs = first_item.wcs  # Assume all elements share the same WCS
    elif isinstance(collection, np.ndarray):
        if collection.ndim == 3 and collection.shape[0] in [1, 3]:  # Grayscale or RGB image
            ax_count = 1
            wcs = None  # No WCS for raw numpy arrays
        elif collection.ndim == 3:  # Multi-channel data (N, H, W)
            ax_count = collection.shape[0]
            wcs = None
        else:
            raise ValueError("Input ndarray must have shape (N, H, W) or (3, H, W) for color images.")
    else:
        raise TypeError("collection must be an NDCollection, a 3D NumPy ndarray, or a color_image array.")

    if not isinstance(vmin, list):
        vmin = [vmin for _ in range(ax_count)]
    if not isinstance(vmax, list):
        vmax = [vmax for _ in range(ax_count)]

    if lat_ticks is None:
        lat_ticks = np.arange(-90, 90, 2) * u.degree
    if lon_ticks is None:
        lon_ticks = np.arange(-180, 180, 2) * u.degree

    fig, axs = plt.subplots(nrows=1, ncols=ax_count, figsize=figsize, sharey=True,
                            subplot_kw={"projection": wcs} if wcs else {})
    if ax_count == 1:
        axs = [axs]

    for i in range(ax_count):
        if isinstance(collection, dict):
            this_cube = collection[collection_keys[i]]
            this_cube.plot(axes=axs[i], cmap=cmap, vmin=vmin[i], vmax=vmax[i])
            im = axs[i].get_images()[0]
            axs[i].set_title(f"{this_cube.meta['POLAR']} at {this_cube.meta['DATE-OBS'][0:16]}")
        elif collection.ndim == 3 and collection.shape[0] == 3:  # RGB image
            im = axs[i].imshow(np.moveaxis(collection, 0, -1), **kwargs)
        else:
            im = axs[i].imshow(collection[i], cmap=cmap, vmin=vmin[i], vmax=vmax[i], **kwargs)
        if wcs:
            axs[i].coords[0].set_ticks(lon_ticks)
            axs[i].coords[1].set_ticks(lat_ticks)
            axs[i].coords[0].set_major_formatter(major_formatter)
            axs[i].coords[1].set_major_formatter(major_formatter)
        axs[i].set_xlabel(xlabel, fontsize=fontsize)
        axs[i].set_ylabel(ylabel, fontsize=fontsize)
        axs[i].tick_params(axis="both", labelsize=fontsize)
        axs[i].grid(color="white", ls="dotted")

    if show_colorbar:
        fig.colorbar(im, orientation="horizontal", ax=axs, shrink=0.9)

    return fig, axs


def get_colormap_str(meta: fits.Header) -> str:
    """Retrieve a color map name from an input FITS file.

    Parameters
    ----------
    meta : fits.Header
        header of the data

    Returns
    -------
    str
        name of appropriate colormap
    """
    if meta["INSTRUME"] == "LASCO":
        detector_name = meta["DETECTOR"]
        if "C2" in detector_name:
            color_map = "soholasco2"
        elif "C3" in detector_name:
            color_map = "soholasco3"
        else:
            warnings.warn("No valid instrument found, setting color_map soholasco2")
            color_map = "soholasco2"
    elif meta["INSTRUME"] == "COSMO K-Coronagraph":
        color_map = "kcor"
    elif meta["INSTRUME"] == "SECCHI":
        detector_name = meta["DETECTOR"]
        if "COR1" in detector_name:
            color_map = "stereocor1"
        elif "COR2" in detector_name:
            color_map = "stereocor2"
        else:
            warnings.warn("No valid instrument found, setting color_map soholasco2")
            color_map = "soholasco2"
    else:
        warnings.warn("No valid instrument found, setting color_map soholasco2")
        color_map = "soholasco2"

    return color_map



def generate_rgb_image(collection,
                        enhancement_method='nrgf',
                        mask_params=None,
                        enhancement_params=None):
    """
    Generate an RGB color image from an NDCollection based on Patel et al. 2023 Res. Notes AAS 7 241.

    Parameters:
    ----------
    collection : NDCollection
        A collection of NDCube objects containing solar data.
    enhancement_method : str
        The radial enhancement method to use. Can be 'intensity_enhance', 'nrgf', 'fnrgf', 'rhef', or 'none'.
        Default is 'nrgf'.
    mask_params : dict, optional
        Dictionary of masking parameters for inner and outer radius. Default values are used if not provided.
        Example:
        - {'inner_radius': 3, 'outer_radius': 32}
    enhancement_params: dict, optional
        Dictionary of parameters specific to enhancement method above.
        Example (check sunkit_image.radial for further parameters information):
        - For 'intensity_enhance': {'scale': 1, 'degree': 1}
        - For 'nrgf': {'inner_radius': 1, 'outer_radius': 32, 'mask_radius': 6}
        - For 'fnrgf': {'order': 3, 'number_angular_segment': 130}
        - For 'rhef': {'supsilon': 0.35, 'fill': np.nan}
        - Empty or None for 'none'.

    Returns:
    -------
    np.ndarray
        Generated color image with RGB channels.
    """
    if mask_params is None:
        mask_params = {}  # Default to an empty dictionary

    if enhancement_params is None:
        enhancement_params = {}  # Default to an empty dictionary

    # Extract mask parameters with defaults
    inner_radius = mask_params.get('inner_radius', 3)
    outer_radius = mask_params.get('outer_radius', 32)

    out_cube = []
    collection_keys = list(collection.keys())

    radial_bin_edges = equally_spaced_bins(inner_radius, outer_radius, collection[collection_keys[0]].data.shape[0] // 4)
    radial_bin_edges *= u.R_sun

    # Define the enhancement function based on the selected method
    enhancement_methods = {
        'intensity_enhance': intensity_enhance,
        'nrgf': nrgf,
        'fnrgf': fnrgf,
        'rhef': rhef,
        'none': lambda x: x  # No enhancement
    }

    enhancement_func = enhancement_methods.get(enhancement_method)

    if enhancement_method not in enhancement_methods:
        raise ValueError("Invalid enhancement method. Choose 'intensity_enhance', 'nrgf', 'fnrgf', 'rhef', or 'none'.")

    for key in collection_keys:
        inputmap = sunpy.map.Map(collection[key].data, collection[key].wcs)

        if enhancement_func:
            # Apply the selected enhancement method
            enhanced = enhancement_func(inputmap, radial_bin_edges=radial_bin_edges, **enhancement_params)
            masked_enhanced = np.ma.array(enhanced.data, mask=np.isnan(enhanced.data))
        else:
            # No enhancement, use the original data
            masked_enhanced = np.ma.array(inputmap.data, mask=np.isnan(inputmap.data))

        scaled = (np.clip(masked_enhanced, 0, 1) * 255).astype('uint8')
        out_cube.append((key, NDCube(data=scaled, meta=collection[key].meta, wcs=collection[key].wcs)))

    outputs = NDCollection(out_cube, meta={}, aligned_axes="all")
    size_im = (scaled.shape[1], scaled.shape[0])
    color_image = np.zeros((3, size_im[1], size_im[0]), dtype=np.uint8)

    color_image[0, :, :] = outputs['Z'].data
    color_image[1, :, :] = outputs['M'].data
    color_image[2, :, :] = outputs['P'].data

    return color_image
