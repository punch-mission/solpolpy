# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest
from pytest import fixture

import astropy.units as u

from solpolpy.core import determine_input_kind, resolve
from solpolpy.polarizers import mzp_to_bpb, bpb_to_btbr


@fixture
def example_mzp():
    d = {-60*u.degree: np.array([2]), 0*u.degree: np.array([1]), 60*u.degree: np.array([1.5]),
         "alpha": np.array([5])*u.degree}
    return d


def test_determine_input_kind_mzp(example_mzp):
    assert determine_input_kind(example_mzp), "MZP"


def test_determine_input_kind():
    d = {"B": np.array([0]), "pB": np.array([1])}
    determined_kind = determine_input_kind(d)
    assert determined_kind, "BpB"


def test_determine_input_kind_fail():
    d = {"B": np.array([0]), "M": np.array([1])}
    with pytest.raises(ValueError):
        determine_input_kind(d)


def test_bpb_to_btbr():
    d = {"B": np.array([2]), "pB": np.array([1])}
    result = resolve(d, "BtBr")
    assert isinstance(result, dict)


def test_mzp_to_btbr(example_mzp):
    result = resolve(example_mzp, "BtBr")
    assert isinstance(result, dict)
