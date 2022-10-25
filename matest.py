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

def plot_BpB_lasco():
    file_list=[path_to_test_files+"lasco_0.fts",
           path_to_test_files+"lasco_-60.fts",
           path_to_test_files+"lasco_+60.fts"]

    output=sp.resolve(file_list, 'BpB')

    minval=np.nanpercentile(output['pB'], 0.0)
    maxval=np.nanpercentile(output['pB'], 100.0)

    figure, ax=plt.subplots()
    ax.imshow(output['pB'], vmin=minval, vmax=maxval)
    plt.show()


def plot_BpB_stereo():
    file_list=[path_to_test_files+"stereo1.fts",
           path_to_test_files+"stereo2.fts",
           path_to_test_files+"stereo3.fts"]


    output=sp.resolve(file_list, 'BpB')

    minval=np.nanpercentile(output['pB'], 5.0)
    maxval=np.nanpercentile(output['pB'], 99.0)

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


def plot_BpB_Gamera():
    file_list=[path_to_test_files+"PBfor_DOICMEM_00000dens_only_gamera_format_0070.fits",
           path_to_test_files+"TBfor_DOICMEM_00000dens_only_gamera_format_0070.fits"]

    data_out={}
    
    hdul = fits.open(file_list[0])
    pBdata_synth=hdul[0].data
    pBdata_synth[pBdata_synth==-9999.0]=np.nan
    data_out['pB']=pBdata_synth

    hdul = fits.open(file_list[1])
    tBdata_synth=hdul[0].data
    tBdata_synth[tBdata_synth==-9999.0]=np.nan
    data_out['B']=tBdata_synth
    ratioData_synth=pBdata_synth/tBdata_synth

    output=sp.resolve(data_out, 'MZP')

    #minval=np.nanpercentile(output_plot, 0.0)
    #maxval=np.nanpercentile(output_plot, 100.0)

    P_synth=output[60*u.degree]
    M_synth=output[-60*u.degree]
    Z_synth=output[0*u.degree]


    figure, ax=plt.subplots(2, 2)
    fontsize=12
    ax[0, 0].imshow(ratioData_synth, vmin=np.nanpercentile(ratioData_synth, 5), vmax=np.nanpercentile(ratioData_synth, 99))
    ax[0, 0].set_title('Ratio', fontsize=fontsize-3)
    ax[0, 1].imshow(P_synth, vmin=np.nanpercentile(P_synth, 10), vmax=np.nanpercentile(P_synth, 90))
    ax[0, 1].set_title('P_synth', fontsize=fontsize-3)
    ax[1, 0].imshow(M_synth, vmin=np.nanpercentile(M_synth, 10), vmax=np.nanpercentile(P_synth, 90))
    ax[1, 0].set_title('M_synth', fontsize=fontsize-3)
    ax[1, 1].imshow(Z_synth, vmin=np.nanpercentile(Z_synth, 10), vmax=np.nanpercentile(P_synth, 90))
    ax[1, 1].set_title('Z_synth', fontsize=fontsize-3)

    plt.show()





plot_MZP()
#plot_BpB_lasco()
#plot_BpB_stereo()
#plot_BpB_Gamera()



#hdul=fits.open(file_list[0])


#minval=np.nanpercentile(hdul[0].data, 0.0)
#maxval=np.nanpercentile(hdul[0].data, 100.0)

#figure, ax=plt.subplots()
#ax.imshow(hdul[0].data, vmin=minval, vmax=maxval)
#plt.show()