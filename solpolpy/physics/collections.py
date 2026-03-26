"""Helpers for building and reading NDCollection-based transform inputs."""

from __future__ import annotations

import copy

from ndcube import NDCollection

from solpolpy.util import combine_all_collection_masks


def data_keys(collection: NDCollection) -> list[str]:
    """Return collection keys excluding the auxiliary alpha cube."""
    return [key for key in collection.keys() if key != "alpha"]


def template_cube(collection: NDCollection, preferred_key: str | None = None):
    """Return the cube that should be used as metadata/WCS template."""
    if preferred_key is not None and preferred_key in collection:
        return collection[preferred_key]
    return collection[data_keys(collection)[0]]


def stack_data(collection: NDCollection, keys: list[str]) -> list:
    """Stack raw cube data in key order."""
    return [collection[key].data for key in keys]


def clone_meta(cube, **updates):
    """Deep copy a cube metadata mapping and update it."""
    meta = copy.deepcopy(cube.meta)
    meta.update(updates)
    return meta


def combine_mask(collection: NDCollection):
    """Combine masks across all data cubes in a collection."""
    return combine_all_collection_masks(collection)
