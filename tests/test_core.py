# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest
from solpolpy import core

#Calculate B_theta using BR, BT



#check output is the same as input for simple array
@pytest.mark.parametrize("inputarray, outputarray",
                         [
                         ({'B':5,'pB':6},{'B':5,'pB':6}),
                         ({'B':1,'pB':2},{'B':1,'pB':2})
                         ])

def test_input_output(inputarray,outputarray):
    np.testing.assert_equal(core.data_constructor(inputarray) ,outputarray)
