# LASCO BpB -> MZP
import numpy as np
import solpolpy as sp
import matplotlib.pyplot as plt
import astropy.units as u

# plot defaults
fontsize=18

# get test data
file_path="../tests/test_support_files/"

lasco_file_list_BpB=[file_path+"lasco_clear.fts",
                    file_path+"lasco_+60.fts"]
# P worked for -60
# Z worked for 0, with a x of weird values
# M worked for +60


output=sp.resolve(lasco_file_list_BpB, 'MZP', alpha='radial90')

M_data=output[-60*u.degree]
M_minval=np.nanpercentile(M_data, 5.0)
M_maxval=np.nanpercentile(M_data, 95.0)

Z_data=output[0*u.degree]
Z_minval=np.nanpercentile(Z_data, 5.0)
Z_maxval=np.nanpercentile(Z_data, 95.0)

P_data=output[60*u.degree]
P_minval=np.nanpercentile(P_data, 5.0)
P_maxval=np.nanpercentile(P_data, 95.0)


figure, ax=plt.subplots(1, 3)
figure.suptitle("LASCO M, Z, P", fontsize=fontsize)

ax[0].imshow(M_data, vmin=M_minval, vmax=M_maxval)
ax[0].set_title('M', fontsize=fontsize-3)

ax[1].imshow(Z_data, vmin=Z_minval, vmax=Z_maxval)
ax[1].set_title('Z', fontsize=fontsize-3)

ax[2].imshow(P_data, vmin=P_minval, vmax=P_maxval)
ax[2].set_title('P', fontsize=fontsize-3)
figure.show()

print()