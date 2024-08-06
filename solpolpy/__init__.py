import importlib

from solpolpy.core import resolve
from solpolpy.instruments import load_data
from solpolpy.plotting import get_colormap_str, plot_collection
from solpolpy.transforms import System

__version__ = importlib.metadata.version("solpolpy")
__all__ = [resolve, load_data, get_colormap_str, plot_collection, System]
