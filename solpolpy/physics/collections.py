"""Helpers for building and reading NDCollection-based transform inputs."""

import copy

import astropy.units as u
from __future__ import annotations
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
    normalized_updates = {}
    for key, value in updates.items():
        if isinstance(value, u.Quantity):
            normalized_updates[key] = value.value
        else:
            normalized_updates[key] = value

    meta = copy.deepcopy(cube.meta)
    if hasattr(meta, "update"):
        meta.update(normalized_updates)
        return meta

    if hasattr(cube.meta, "to_fits_header"):
        meta = cube.meta.to_fits_header(cube.wcs)
        meta.update(normalized_updates)
        return meta

    meta = {key: copy.deepcopy(cube.meta[key]) for key in cube.meta.keys()}
    meta.update(normalized_updates)
    return meta


def combine_mask(collection: NDCollection):
    """Combine masks across all data cubes in a collection."""
    return combine_all_collection_masks(collection)
