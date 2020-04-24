import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from astropy.io import fits
from astropy import units
from astropy import constants as const
import scipy.integrate as integrate
from matplotlib import animation
from numpy.random import randint
import matplotlib.gridspec as gridspec
from scipy import ndimage

plt.ion()

cm_inferno = plt.get_cmap('inferno')
cm_inferno.set_bad(color='black')
imkw_inferno = dict(origin='lower', interpolation='none', cmap=cm_inferno)

C = 1.82243e18 

#SIMU SAURY THIN MAMD
path = "/data/amarchal/ROHSA_paper/data/synthetic_obs/"
filename = "Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA.fits"
hdu = fits.open(path + filename) 
hdr = hdu[0].header
cube = hdu[0].data
cube[np.where(cube != cube)] = 0.

reso = 0.8
noise = 0.05

# Smooth with 3D gaussian filter
smooth_cube_obs = ndimage.filters.gaussian_filter(cube, [0,2,2])

# #Add noise in the synthetic cube
for i in range(cube.shape[1]):
    for j in range(cube.shape[2]):
        smooth_cube_obs[:, i, j] += np.random.randn(smooth_cube_obs.shape[0]) * noise

path_out = '/data/amarchal/ROHSA_paper/data/synthetic_obs/'
fileout = filename[:-5] + '_noise_' + str(noise) + '_K' + '_beam_0_2_2.fits'

# Write the observation cube
hdu0 = fits.PrimaryHDU(smooth_cube_obs)
# hdu0.header['COMMENT'] = 'Brightness Temperature Tb / noise '  + str(noise) + ' K'
hdu0.header['COMMENT'] = 'Brightness Temperature Tb'
hdu0.header['NAXIS'] = 3
hdu0.header['NAXIS1'] = 256
hdu0.header['NAXIS2'] = 256
hdu0.header['NAXIS3'] = cube.shape[0]
hdu0.header['CTYPE3'] = 'v [km.s-1]'
hdu0.header['CRVAL3'] = hdr['CRVAL3']
hdu0.header['CDELT3'] = hdr['CDELT3']
hdu0.header['CRPIX3'] = hdr['CRPIX3']
hdu0.header['BUNIT'] = 'K'

hdulist = fits.HDUList([hdu0])
hdulist.writeto(path_out + fileout, overwrite=True) 

