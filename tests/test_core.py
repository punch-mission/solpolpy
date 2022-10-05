# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest

from solpolpy.core import determine_input_kind, resolve
from solpolpy.polarizers import mzp_to_bpb, bpb_to_btbr

def test_true():
    assert True


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


def test_mzp_to_btbr():
    d = {np.deg2rad(-60): np.array([2]), np.deg2rad(0): np.array([1]), np.deg2rad(60): np.array([1.5]),
         "alpha": np.array([np.deg2rad(5)])}
    result = resolve(d, "BtBr")
    assert isinstance(result, dict)
