import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from astropy.io import fits
from astropy import units
from astropy import constants as const
from astropy import wcs

# from mylibrary.CGD import mod_tools
# from mylibrary.PROJ import mod_projection
# from mylibrary.PLOT import mod_plot

plt.ion()
cm = plt.get_cmap('inferno')
cm.set_bad(color='black')
imkw = dict(origin='lower', interpolation='none', cmap=cm)

# Constant
m_h = 1.6737236e-27 #kg
C = 1.83e18 #K-1cm-2 / (km.s-1)
pc2cm = units.pc.to(units.m) * 1.e2

path = '/data/glx-calcul3/data1/amarchal/IRAM/to_decompose/data/'
filename_data = 'Draco9_CO10.fits'

hdu_list_data = fits.open(path + filename_data)
hdu_data = hdu_list_data[0]

if len(hdu_list_data[0].shape) == 4.:
        cube = hdu_data.data[0]
elif len(hdu_list_data[0].shape) == 3 :
        cube = hdu_data.data
else: print 'Invalid dimension of input data'

cube[np.where(cube != cube)] = 0

subcube = cube

nx = subcube.shape[2]
ny = subcube.shape[1]
nv = subcube.shape[0]

with open(path + filename_data[:-5] + '.dat','w+') as f:
        f.write('{:d}\t{:d}\t{:d}\n'.format(nv, ny, nx))
        for i in range(ny):
                for j in range(nx):
                        for k in range(nv):
                                line = '{:d}\t{:d}\t{:d}\t{:0.16f}\n'.format(k, i, j, subcube[k,i,j])
                                f.write(line)
                
