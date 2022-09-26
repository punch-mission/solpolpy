# -*- coding: utf-8 -*-
"""
pytest test suite for the polarizers module of solpolpy
"""

import numpy as np
import pytest
from solpolpy import validation as valid


#check output is the same as input for simple array
@pytest.mark.parametrize("inputarray, outputarray",
                         [(np.array([[1, 2, 3], 
                                     [4, 5, 6], 
                                     [7, 8, 9]]),
                           np.array([[1, 2, 3], 
                                     [4, 5, 6], 
                                     [7, 8, 9]]))])

def test_input_output(inputarray,outputarray):
    np.testing.assert_array_equal(valid.data_type(inputarray) ,[inputarray])

