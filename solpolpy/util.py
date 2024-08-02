import numpy as np
from ndcube import NDCollection

from solpolpy.errors import MissingAlphaError


def combine_all_collection_masks(collection: NDCollection) -> np.ndarray | None:
    """Combine all the masks in a given collection."""
    return combine_masks(*[cube.mask for key, cube in collection.items() if key != "alpha"])


def combine_masks(*args) -> np.ndarray | None:
    """Combine masks.

    If any of the masks are None, the result is None.
    Otherwise, when combining any value that is masked in any of the input args, gets masked, i.e. it does a logical or.
    """
    if any(arg is None for arg in args):
        return None
    else:
        return np.logical_or.reduce(args)


def use_alpha(transform):
    """Require a transform to use the alpha parameter and check for its existence."""

    def wrapped(input_collection, *func_args, **func_kwargs):
        if "alpha" not in input_collection:
            msg = "alpha expected in input_collection but not found."
            raise MissingAlphaError(msg)
        return transform(input_collection, *func_args, **func_kwargs)
    wrapped.uses_alpha = True

    return wrapped
