import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
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

def plot_spect(x,y,cube,id1,id2,amp,mean,sigma,v):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.step(v, cube[:,x,y], color='cornflowerblue', linewidth=1.)
    ax.plot(v, reconstructed_cube[:,x,y], color='r')
    iid = np.where((id1 == x) & (id2 == y))[0]
    for i in range(len(iid)):
        ax.plot(v, gaussian(np.arange(cube.shape[0]), amp[iid[i]], mean[iid[i]], sigma[iid[i]]))
    for i in range(len(iid)):
        print amp[iid[i]], mean[iid[i]], sigma[iid[i]]

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

#FULL FIELD
data = np.genfromtxt("/data/amarchal/ROHSA_paper/ROHSA/GHIGLS_NEP_Tb_gauss_run_4.dat") # 4

id1 = data[:, 0]
id2 = data[:, 1]
sigma = data[:, 4]
amp = data[:, 2]
mean = data[:, 3] - 1

vmean = mean2vel(CRVAL, CDELT, CRPIX, mean)

reconstructed_cube = np.zeros((dim_v, dim_y, dim_x))

for i in range(len(amp)):
    gauss = gaussian(np.arange(dim_v), amp[i], mean[i], sigma[i])
    reconstructed_cube[:, int(id1[i]), int(id2[i])] += gauss

n_gauss = len(id1) / (dim_y*dim_x)
params = np.zeros((3*n_gauss, dim_y, dim_x))

i__ = 0
for i in range(dim_y):
    for j in range(dim_x):
        for k in range(n_gauss):
            params[0+(3*k),i,j] = amp[i__]
            params[1+(3*k),i,j] = vmean[i__]
            params[2+(3*k),i,j] = sigma[i__] * reso
            i__ += 1
            
fields = [np.sqrt(2.*np.pi) * params[0+(3*k)] * params[2+(3*k)] for k in np.arange(n_gauss)]
field = [f * C / 1.e18 for f in fields]
ampfield = [params[0+(3*k)] for k in np.arange(n_gauss)]
vfield = [params[1+(3*k)] for k in np.arange(n_gauss)]
sigfield = [params[2+(3*k)] for k in np.arange(n_gauss)]

iddx = np.argsort(np.mean(vfield, axis=(1,2)))

field = [field[idd] for idd in iddx]
vfield = [vfield[idd] for idd in iddx]
ampfield = [ampfield[idd] for idd in iddx]
sigfield = [sigfield[idd] for idd in iddx]

ys = 100
# xs = np.arange(200)+10
xs = np.arange(198)

#Plot mosaic spectra random walk 2D 
pvalues = np.logspace(-0.7, 0, 12)
pmin = pvalues[0]
pmax = pvalues[-1]

def norm(pval):
    return (pval - pmin) / float(pmax - pmin)

for l in xs:      
    ny = 4; nx = 4                                                                                                                                         
    center_y = ys; center_x = xs[l]
    cb = "magenta"
    cw = "crimson"
    fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(cm2inch((18.,14.))))
    fig.subplots_adjust(hspace=0, wspace=0, left=0, right=1, top=1, bottom=0)
    for i in np.arange(ny):
        for j in np.arange(nx):
            # axs[i][j].set_yticks(np.arange(0., 18., 5.))
            # axs[i][j].set_xticks(np.arange(-10, 20., 10.))
            # axs[i][j].set_ylim([-1,20])
            axs[i][j].set_xlim([-80,60])
            axs[i][j].set_ylim([0,12])
            axs[i][j].step(v, cube[:,center_y+i,center_x+j], color='cornflowerblue', linewidth=2.)
            axs[i][j].plot(v, reconstructed_cube[:,center_y+i,center_x+j], linestyle="-", linewidth=2., color="k")
            iid = np.where((id1 == center_y+i) & (id2 == center_x+j))[0]
            for k in range(len(iid)):            
                axs[i][j].plot(v, gaussian(np.arange(len(v)), amp[iid[k]], mean[iid[k]], sigma[iid[k]]), linewidth=2., color=plt.cm.inferno(pvalues[k]))
                # axs[i][j].plot(v, reconstructed_cube_cnm[:,center_y+i,center_x+j], linewidth=1., color=cb)
                # axs[i][j].plot(v, reconstructed_cube_wnm[:,center_y+i,center_x+j], linewidth=1., color=cw)
            if j == 0: axs[i][j].set_ylabel(r'T [k]')
            axs[i][j].set_xlabel(r'v [km s$^{-1}$]')
    plt.savefig("plot/mosaic_spectra/" + 'mosaic_spectra_all_NEP_{}.png'.format(l), format='png', bbox_inches='tight', pad_inches=0.02)

