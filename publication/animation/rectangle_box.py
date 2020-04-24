import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from astropy.io import fits
from astropy import units
from astropy import constants as const
import scipy.integrate as integrate
from matplotlib import animation
from numpy.random import randint
import matplotlib.gridspec as gridspec
from astropy import wcs
from scipy.stats import moment
from scipy import stats
import pandas as pd
import seaborn as sns
from PIL import Image
import pylab 
import random
from matplotlib.patches import Rectangle

import turbulence as tb

# plt.ion()

cm_coolwarm = plt.get_cmap('coolwarm')
cm_coolwarm.set_bad(color='black')
imkw_coolwarm = dict(origin='lower', interpolation='none', cmap=cm_coolwarm)

cm_inferno = plt.get_cmap('inferno')
cm_inferno.set_bad(color='black')
imkw_inferno = dict(origin='lower', interpolation='none', cmap=cm_inferno)

cm_viridis = plt.get_cmap('viridis')
cm_viridis.set_bad(color='black')
imkw_viridis = dict(origin='lower', interpolation='none', cmap=cm_viridis)

cm_cubehelix = plt.get_cmap('cubehelix')
cm_cubehelix.set_bad(color='black')
imkw_cubehelix = dict(origin='lower', interpolation='none', cmap=cm_cubehelix)

def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu)**2)/(2. * sigma**2))

def mean2vel(CRVAL, CDELT, CRPIX, mean):
    vel = [(CRVAL + CDELT * (mean[i] - CRPIX)) for i in range(len(mean))] #VLSR [km.s-1] #FIXME
    return vel

def vel2mean(CRVAL, CDELT, CRPIX, vel):
    return [((vel[i] - CRVAL) / CDELT + CRPIX) for i in range(len(vel))]

def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i/inch for i in tupl[0])
        else:
            return tuple(i/inch for i in tupl)

C = 1.82243e18 

#DATA
hdu = fits.open("/data/amarchal/GHIGLS/data/GHIGLS_NEP_Tb.fits") 
hdr = hdu[0].header

cube = hdu[0].data[0]
cube[np.where(cube != cube)] = 0.

cube = cube[200:500,:,:]

dim_x = cube.shape[2]
dim_y = cube.shape[1]
dim_v = cube.shape[0]

CDELT = hdr['CDELT3'] *1.e-3  #km.s-1                                                                                                      
CRVAL = hdr['CRVAL3'] *1.e-3  #km.s-1                                                                                                             
CRPIX = hdr['CRPIX3'] - 200

reso = np.abs(CDELT)

v = mean2vel(CRVAL, CDELT, CRPIX, np.arange(cube.shape[0]))
original_map = np.sum(cube, axis=0) * reso * C / 1.e18

w = wcs.WCS(naxis=2)
w.wcs.crpix = [hdr['CRPIX1'], hdr['CRPIX2']]
w.wcs.cdelt = np.array([hdr['CDELT1'], hdr['CDELT2']])
w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
w.wcs.ctype = [hdr['CTYPE1'], hdr['CTYPE2']]

for k in np.arange(198):
    rect = Rectangle((6+k,100),6,6,linewidth=1,edgecolor='m',facecolor='cyan', alpha=0.4)
    fig = plt.figure(figsize=[6,6])
    ax = fig.add_subplot(111)
    ax.imshow(original_map, origin="lower", cmap="inferno")
    ax.add_patch(rect)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_frame_on(False)
    plt.savefig("plot/rectangle_box/" + 'NHI_TOT_NEP_rect_{}.png'.format(k), format='png', bbox_inches='tight', pad_inches=0.02)

