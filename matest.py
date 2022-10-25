#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test series of images
"""

import numpy as np
import solpolpy as sp
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import astropy.units as u


TESTDATA_DIR = os.path.dirname(__file__)
path_to_test_files=TESTDATA_DIR+'/tests/test_support_files/'

def plot_BpB():
    file_list=[path_to_test_files+"stereo1.fts",
           path_to_test_files+"stereo2.fts",
           path_to_test_files+"stereo3.fts"]

    file_list=[path_to_test_files+"lasco_0.fts",
           path_to_test_files+"lasco_-60.fts",
           path_to_test_files+"lasco_+60.fts"]

    output=sp.resolve(file_list, 'BpB')

    stereominval=np.nanpercentile(output['pB'], 5.0)
    stereomaxval=np.nanpercentile(output['pB'], 99.0)

    minval=np.nanpercentile(output['pB'], 0.0)
    maxval=np.nanpercentile(output['pB'], 100.0)

    figure, ax=plt.subplots()
    ax.imshow(output['pB'], vmin=minval, vmax=maxval)
    plt.show()


def plot_MZP():
    file_list=[path_to_test_files+"lasco_clear.fts",
           path_to_test_files+"lasco_-60.fts"]
    output=sp.resolve(file_list, 'MZP')

    output_plot=output[-60*u.degree]

    minval=np.nanpercentile(output_plot, 0.0)
    maxval=np.nanpercentile(output_plot, 100.0)

    figure, ax=plt.subplots()
    ax.imshow(output_plot, vmin=minval, vmax=maxval)
    plt.show()


plot_BpB()





#hdul=fits.open(file_list[0])


#minval=np.nanpercentile(hdul[0].data, 0.0)
#maxval=np.nanpercentile(hdul[0].data, 100.0)

#figure, ax=plt.subplots()
#ax.imshow(hdul[0].data, vmin=minval, vmax=maxval)
#plt.show()