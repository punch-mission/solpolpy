import numpy as np
import solpolpy as sp
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import astropy.units as u

# plot defaults
fontsize=18

# get test data
file_path="../tests/test_support_files/"
gamera_file_list_BpB=[file_path+"TBfor_DOICMEM_00000dens_only_gamera_format_0070.fits",
                      file_path+"PBfor_DOICMEM_00000dens_only_gamera_format_0070.fits"]

lasco_file_list_MZP=[file_path+"lasco_+60.fts",
                     file_path+"lasco_-60.fts",
                     file_path+"lasco_0.fts"]

lasco_file_list_BpB=[file_path+"lasco_clear.fts",
                    file_path+"lasco_-60.fts"]

stereo_file_list_MZP=[file_path+"stereo_0.fts",
                      file_path+"stereo_120.fts",
                      file_path+"stereo_240.fts"]

# Gamera BpB -> MZP

# output=sp.resolve(gamera_file_list_BpB, 'MZP')

data_out = {}

hdul = fits.open(gamera_file_list_BpB[1])
pBdata_synth = hdul[0].data
pBdata_synth[pBdata_synth == -9999.0] = np.nan
data_out['pB'] = pBdata_synth

hdul = fits.open(gamera_file_list_BpB[0])
tBdata_synth = hdul[0].data
tBdata_synth[tBdata_synth == -9999.0] = np.nan
data_out['B'] = tBdata_synth
# ratioData_synth=pBdata_synth/tBdata_synth

output_gamera = sp.resolve(data_out, 'MZP', alpha='radial90')
# output_plot=output[-60*u.degree]

M_minval = np.nanpercentile(output_gamera[-60 * u.degree], 5.0)
M_maxval = np.nanpercentile(output_gamera[-60 * u.degree], 99.0)

Z_minval = np.nanpercentile(output_gamera[0 * u.degree], 5.0)
Z_maxval = np.nanpercentile(output_gamera[0 * u.degree], 99.0)

P_minval = np.nanpercentile(output_gamera[60 * u.degree], 5.0)
P_maxval = np.nanpercentile(output_gamera[60 * u.degree], 99.0)

figure, ax = plt.subplots(1, 3)
figure.suptitle("GAMERA M, Z, P", fontsize=fontsize)

ax[0].imshow(output_gamera[-60 * u.degree], vmin=M_minval, vmax=M_maxval)
ax[0].set_title('M', fontsize=fontsize - 3)

ax[1].imshow(output_gamera[0 * u.degree], vmin=Z_minval, vmax=Z_maxval)
ax[1].set_title('Z', fontsize=fontsize - 3)

ax[2].imshow(output_gamera[60 * u.degree], vmin=P_minval, vmax=P_maxval)
ax[2].set_title('P', fontsize=fontsize - 3)
figure.show()
