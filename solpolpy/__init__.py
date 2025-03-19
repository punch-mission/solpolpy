import importlib

from solpolpy.core import resolve
from solpolpy.instruments import load_data
from solpolpy.plotting import generate_rgb_image, get_colormap_str, plot_collection
from solpolpy.transforms import System
from solpolpy.util import collection_to_maps

__version__ = importlib.metadata.version("solpolpy")
__all__ = [resolve, load_data, get_colormap_str, plot_collection, generate_rgb_image, collection_to_maps, System]
