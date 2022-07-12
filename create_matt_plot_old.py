
"""Docstring for the math_equations.py module.

Created by Bryce M. Walbridge on 6/02/2022.

The module explored definitions of math equation illistrated by the article, 
Three-polarizer Treatment of Linear Polarization in Coronagraphs and 
Heliospheric Imagers by Craig E. DeForester, Daniel B. Seaton, and 
Matthew J. West. The overall goal is to convert values from three polarized 
images (B, pB pB') or Stokes (I, Q, U) representations of linear
polarization from polarizer triplet data. This based on the formulae proposed
in this article.

"""

import numpy as np 
import matplotlib.pyplot as plt

# create plotting functions

def create_plot(data, title='title', lower_percentile=5, upper_percentile=99.0, fontsize=16, nan=None):
    '''
    Create a plot between lower and upper percentile for ordinary data
    '''
    if nan != None:
        plotmin=np.nanpercentile(data, lower_percentile)
        plotmax=np.nanpercentile(data, upper_percentile)
    else:
        plotmin=np.percentile(data, lower_percentile)
        plotmax=np.percentile(data, upper_percentile)

    plt.figure()
    plt.imshow(data, vmin=plotmin, vmax=plotmax)
    plt.title(title, fontsize=fontsize)
    plt.tight_layout()
    plt.show()


def create_2_plot(data1, data2, suptitle='big_title', fontsize=16,
                  title1='title1', title2='title2',
                  lower1=5, upper1=99.0, lower2=None, upper2=None,
                  nan1=None, nan2=None,
                  abs1=None, abs2=None):
    '''
    Create 4 plots between lower and upper percentile for ordinary data
    '''

    fig, axs = plt.subplots(1, 2)
    fig.suptitle(suptitle, fontsize=fontsize)

    if abs1 != None:
        plotmin1=lower1
        plotmax1=upper1
    elif nan1 != None:
        plotmin1=np.nanpercentile(data1, lower1)
        plotmax1=np.nanpercentile(data1, upper1)
    else:
        plotmin1=np.percentile(data1, lower1)
        plotmax1=np.percentile(data1, upper1)

    axs[0].imshow(data1, vmin=plotmin1, vmax=plotmax1)
    axs[0].set_title(title1, fontsize=fontsize-2)

    if lower2==None:
        lower2=lower1 
    if upper2==None:
        upper2=upper1
    if abs2 != None:
        plotmin2=lower2
        plotmax2=upper2
    elif nan2 != None:
        plotmin2=np.nanpercentile(data2, lower2)
        plotmax2=np.nanpercentile(data2, upper2)
    else:
        plotmin2=np.percentile(data2, lower2)
        plotmax2=np.percentile(data2, upper2)

    axs[1].imshow(data2, vmin=plotmin2, vmax=plotmax2)
    axs[1].set_title(title2, fontsize=fontsize-2)

    fig.tight_layout()
    fig.show()


def create_3_plot(data1, data2, data3, suptitle='big_title', fontsize=16,
                  title1='title1', title2='title2', title3='title3',
                  lower1=5, upper1=99.0, lower2=None, upper2=None,
                  lower3=None, upper3=None,
                  nan1=None, nan2=None, nan3=None,
                  abs1=None, abs2=None, abs3=None):
    '''
    Create 4 plots between lower and upper percentile for ordinary data
    '''

    fig, axs = plt.subplots(1, 3)
    fig.suptitle(suptitle, fontsize=fontsize)

    if abs1 != None:
        plotmin1=lower1
        plotmax1=upper1
    elif nan1 != None:
        plotmin1=np.nanpercentile(data1, lower1)
        plotmax1=np.nanpercentile(data1, upper1)
    else:
        plotmin1=np.percentile(data1, lower1)
        plotmax1=np.percentile(data1, upper1)

    axs[0].imshow(data1, vmin=plotmin1, vmax=plotmax1)
    axs[0].set_title(title1, fontsize=fontsize-2)

    if lower2==None:
        lower2=lower1 
    if upper2==None:
        upper2=upper1
    if abs2 != None:
        plotmin2=lower2
        plotmax2=upper2
    elif nan2 != None:
        plotmin2=np.nanpercentile(data2, lower2)
        plotmax2=np.nanpercentile(data2, upper2)
    else:
        plotmin2=np.percentile(data2, lower2)
        plotmax2=np.percentile(data2, upper2)

    axs[1].imshow(data2, vmin=plotmin2, vmax=plotmax2)
    axs[1].set_title(title2, fontsize=fontsize-2)

    if lower3==None:
        lower3=lower1 
    if upper3==None:
        upper3=upper1

    if abs3 != None:
        plotmin3=lower3
        plotmax3=upper3
    elif nan3 != None:
        plotmin3=np.nanpercentile(data3, lower3)
        plotmax3=np.nanpercentile(data3, upper3)
    else:
        plotmin3=np.percentile(data3, lower3)
        plotmax3=np.percentile(data3, upper3)

    axs[2].imshow(data3, vmin=plotmin3, vmax=plotmax3)
    axs[2].set_title(title3, fontsize=fontsize-2)


    fig.tight_layout()
    fig.show()



    

def create_4_plot(data1, data2, data3, data4, suptitle='big_title', fontsize=16,
                  title1='title1', title2='title2', title3='title3', title4='title4',
                  lower1=5, upper1=99.0, lower2=None, upper2=None,
                  lower3=None, upper3=None, lower4=None, upper4=None,
                  nan1=None, nan2=None, nan3=None, nan4=None,
                  abs1=None, abs2=None, abs3=None, abs4=None):
    '''
    Create 4 plots between lower and upper percentile for ordinary data
    '''

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(suptitle, fontsize=fontsize)

    if abs1 != None:
        plotmin1=lower1
        plotmax1=upper1
    elif nan1 != None:
        plotmin1=np.nanpercentile(data1, lower1)
        plotmax1=np.nanpercentile(data1, upper1)
    else:
        plotmin1=np.percentile(data1, lower1)
        plotmax1=np.percentile(data1, upper1)

    axs[0, 0].imshow(data1, vmin=plotmin1, vmax=plotmax1)
    axs[0, 0].set_title(title1, fontsize=fontsize-3)

    if lower2==None:
        lower2=lower1 
    if upper2==None:
        upper2=upper1
    if abs2 != None:
        plotmin2=lower2
        plotmax2=upper2
    elif nan2 != None:
        plotmin2=np.nanpercentile(data2, lower2)
        plotmax2=np.nanpercentile(data2, upper2)
    else:
        plotmin2=np.percentile(data2, lower2)
        plotmax2=np.percentile(data2, upper2)

    axs[0, 1].imshow(data2, vmin=plotmin2, vmax=plotmax2)
    axs[0, 1].set_title(title2, fontsize=fontsize-3)

    if lower3==None:
        lower3=lower1 
    if upper3==None:
        upper3=upper1

    if abs3 != None:
        plotmin3=lower3
        plotmax3=upper3
    elif nan3 != None:
        plotmin3=np.nanpercentile(data3, lower3)
        plotmax3=np.nanpercentile(data3, upper3)
    else:
        plotmin3=np.percentile(data3, lower3)
        plotmax3=np.percentile(data3, upper3)

    axs[1, 0].imshow(data3, vmin=plotmin3, vmax=plotmax3)
    axs[1, 0].set_title(title3, fontsize=fontsize-3)

    if lower4==None:
        lower4=lower1 
    if upper4==None:
        upper4=upper1
    if abs4 != None:
        plotmin4=lower4
        plotmax4=upper4
    elif nan4 != None:
        plotmin4=np.nanpercentile(data4, lower4)
        plotmax4=np.nanpercentile(data4, upper4)
    else:
        plotmin4=np.percentile(data4, lower4)
        plotmax4=np.percentile(data4, upper4)

    axs[1, 1].imshow(data4, vmin=plotmin4, vmax=plotmax4)
    axs[1, 1].set_title(title4, fontsize=fontsize-3)

    fig.tight_layout()
    fig.show()



def create_nan_plot(data, title='title', lower_percentile=5, upper_percentile=99.0, fontsize=16):
    '''
    Create a plot between lower and upper percentile for data with nan's
    '''
    plotmin=np.nanpercentile(data, lower_percentile)
    plotmax=np.nanpercentile(data, upper_percentile)

    plt.figure()
    plt.imshow(data, vmin=plotmin, vmax=plotmax)
    plt.title(title, fontsize=fontsize)
    plt.tight_layout()
    plt.show()


def create_4_plot_nan(data1, data2, data3, data4, suptitle='big_title',
                  title1='title1', title2='title2', title3='title3', title4='title4',
                  lower_percentile=5, upper_percentile=99.0, fontsize=16):
    '''
    Create 4 plots between lower and upper percentile for ordinary data
    '''
        
    fig, axs = plt.subplots(2, 2)
    fig.suptitle(suptitle, fontsize=fontsize)
    
    plotmin=np.nanpercentile(data1, lower_percentile)
    plotmax=np.nanpercentile(data1, upper_percentile)
    axs[0, 0].imshow(data1, vmin=plotmin, vmax=plotmax)
    axs[0, 0].set_title(title1, fontsize=fontsize)

    plotmin=np.nanpercentile(data2, lower_percentile)
    plotmax=np.nanpercentile(data2, upper_percentile)
    axs[0, 1].imshow(data2, vmin=plotmin, vmax=plotmax)
    axs[0, 1].set_title(title2, fontsize=fontsize)

    plotmin=np.nanpercentile(data3, lower_percentile)
    plotmax=np.nanpercentile(data3, upper_percentile)
    axs[1, 0].imshow(data3, vmin=plotmin, vmax=plotmax)
    axs[1, 0].set_title(title3, fontsize=fontsize)

    plotmin=np.nanpercentile(data4, lower_percentile)
    plotmax=np.nanpercentile(data4, upper_percentile)
    axs[1, 1].imshow(data4)
    axs[1, 1].set_title(title4, fontsize=fontsize)
    fig.tight_layout()
    fig.show()