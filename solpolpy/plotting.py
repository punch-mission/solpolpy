import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.visualization.colormaps as cm  # noqa: F401
from astropy.io import fits
from ndcube import NDCollection


def plot_collection(collection: NDCollection,
                    figsize=(8, 8),
                    show_colorbar=False,
                    lat_ticks=None,
                    lon_ticks=None,
                    major_formatter='dd',
                    xlabel="HP Longitude",
                    ylabel="HP Latitude",
                    vmin=None,
                    vmax=None,
                    cmap='Greys_r',
                    ignore_alpha=True,
                    fontsize=18):
    """Plot a solpolpy NDCollection input or output

    Parameters
    ----------
    collection : NDCollection
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

    Returns
    -------
    Matplotlib figure and axes
        the plotted figure and axes are returned for any additional edits
    """
    ax_count = len(collection)
    collection_keys = list(collection.keys())

    if ignore_alpha:
        collection_keys = [k for k in collection_keys if k != "alpha"]
        ax_count = len(collection_keys)

    if not isinstance(vmin, list):
        vmin = [vmin for _ in range(ax_count)]

    if not isinstance(vmax, list):
        vmax = [vmax for _ in range(ax_count)]

    if lat_ticks is None:
        lat_ticks = np.arange(-90, 90, 2) * u.degree

    if lon_ticks is None:
        lon_ticks = np.arange(-180, 180, 2) * u.degree

    fig, axs = plt.subplots(nrows=1,
                            ncols=ax_count,
                            figsize=figsize,
                            sharey=True,
                            subplot_kw={'projection': collection[collection_keys[0]].wcs})
    for i in range(ax_count):
        this_cube = collection[collection_keys[i]]
        this_cube.plot(axes=axs[i], cmap=cmap, vmin=vmin[i], vmax=vmax[i])
        im = axs[i].get_images()[0]
        axs[i].set_title(f"{this_cube.meta['POLAR']} at {this_cube.meta['DATE-OBS'][0:16]}")
        axs[i].coords[0].set_ticks(lon_ticks)
        axs[i].coords[1].set_ticks(lat_ticks)
        axs[i].coords[0].set_major_formatter(major_formatter)
        axs[i].coords[1].set_major_formatter(major_formatter)
        axs[i].set_xlabel(xlabel, fontsize=fontsize)
        axs[i].set_ylabel(ylabel, fontsize=fontsize)
        axs[i].tick_params(axis='both', labelsize=fontsize)
        axs[i].grid(color='white', ls='dotted')

    if show_colorbar:
        fig.colorbar(im, orientation='horizontal', ax=axs, shrink=0.9)

    return fig, axs


def get_colormap_str(meta: fits.Header) -> str:
    """Retrieve a color map name from an input FITS file
    Parameters
    ----------
    meta : fits.Header
        header of the data

    Returns
    -------
    str
        name of appropriate colormap
    """

    if meta['INSTRUME'] == 'LASCO':
        detector_name = meta['DETECTOR']
        if 'C2' in detector_name:
            color_map = 'soholasco2'
        elif 'C3' in detector_name:
            color_map = 'soholasco3'
        else:
            warnings.warn("No valid instrument found, setting color_map soholasco2")
            color_map = 'soholasco2'
    elif meta['INSTRUME'] == 'COSMO K-Coronagraph':
        color_map = 'kcor'
    elif meta['INSTRUME'] == 'SECCHI':
        detector_name = meta['DETECTOR']
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
