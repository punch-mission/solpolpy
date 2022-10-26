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

# LASCO MZP ->BpB

output=sp.resolve(lasco_file_list_MZP, 'BpB')

B_minval=np.nanpercentile(output['B'], 5.0)
B_maxval=np.nanpercentile(output['B'], 99.0)

pB_minval=np.nanpercentile(output['pB'], 5.0)
pB_maxval=np.nanpercentile(output['pB'], 99.0)

figure, ax=plt.subplots(1,2)
figure.suptitle("LASCO B, pB", fontsize=fontsize)

ax[0].imshow(output['B'], vmin=B_minval, vmax=B_maxval)
ax[0].set_title('B', fontsize=fontsize-3)

ax[1].imshow(output['pB'], vmin=pB_minval, vmax=pB_maxval)
ax[1].set_title('pB', fontsize=fontsize-3)

plt.show()

second_trip = sp.resolve(output, "MZP")
figure, ax=plt.subplots(1,3)
figure.suptitle("LASCO M, Z, P", fontsize=fontsize)

ax[0].imshow(second_trip[-60*u.degree])
ax[0].set_title('M', fontsize=fontsize-3)

ax[1].imshow(second_trip[0*u.degree])
ax[1].set_title('Z', fontsize=fontsize-3)

ax[2].imshow(second_trip[60*u.degree])
ax[2].set_title('P', fontsize=fontsize-3)
figure.show()


with fits.open(file_path+"lasco_0.fts") as hdul:
    zero_data = hdul[0].data
print(np.allclose(second_trip[0*u.degree], zero_data))

print()