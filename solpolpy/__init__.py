import importlib

from solpolpy.core import resolve
from solpolpy.instruments import load_data
from solpolpy.plotting import generate_rgb_image, get_colormap_str, plot_collection
from solpolpy.transforms import System
from solpolpy.util import collection_to_maps, solnorth_from_wcs

__version__ = importlib.metadata.version("solpolpy")
__all__ = [
           "System",
           "collection_to_maps",
           "generate_rgb_image",
           "get_colormap_str",
           "load_data",
           "plot_collection",
           "resolve",
           "solnorth_from_wcs",
]
