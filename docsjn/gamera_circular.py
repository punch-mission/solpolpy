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

gamera_in = {}

hdul = fits.open(gamera_file_list_BpB[1])
pBdata_synth = hdul[0].data
pBdata_synth[pBdata_synth == -9999.0] = np.nan
gamera_in['pB'] = pBdata_synth

hdul = fits.open(gamera_file_list_BpB[0])
tBdata_synth = hdul[0].data
tBdata_synth[tBdata_synth == -9999.0] = np.nan
gamera_in['B'] = tBdata_synth

alpha = np.zeros_like(gamera_in['pB'])
output=sp.resolve(gamera_in, 'MZP', alpha="zeros")

figure, ax=plt.subplots(1,3)
figure.suptitle("Gamera M, Z, P", fontsize=fontsize)

ax[0].imshow(output[-60*u.degree])
ax[0].set_title('M', fontsize=fontsize-3)

ax[1].imshow(output[0*u.degree])
ax[1].set_title('Z', fontsize=fontsize-3)

ax[2].imshow(output[60*u.degree])
ax[2].set_title('P', fontsize=fontsize-3)
figure.show()

del output['alpha']
second_trip = sp.resolve(output, "BpB", alpha='radial90')
B_minval=np.nanpercentile(second_trip['B'], 5.0)
B_maxval=np.nanpercentile(second_trip['B'], 99.0)

pB_minval=np.nanpercentile(second_trip['pB'], 5.0)
pB_maxval=np.nanpercentile(second_trip['pB'], 99.0)

figure, ax=plt.subplots(1,2)
figure.suptitle("Gamera B, pB", fontsize=fontsize)

ax[0].imshow(second_trip['B'], vmin=B_minval, vmax=B_maxval)
ax[0].set_title('B', fontsize=fontsize-3)

ax[1].imshow(second_trip['pB'], vmin=pB_minval, vmax=pB_maxval)
ax[1].set_title('pB', fontsize=fontsize-3)

plt.show()


with fits.open(file_path+"PBfor_DOICMEM_00000dens_only_gamera_format_0070.fits") as hdul:
    pb_data = hdul[0].data
print(np.allclose(second_trip["pB"], pb_data))

print()
