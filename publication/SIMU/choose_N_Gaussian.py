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
from astropy import wcs
from scipy.stats import moment
import colorcet as cc

import turbulence as tb 

plt.ion()

cm_coolwarm = cc.cm["coolwarm"]
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

#SIMU SAURY
hdu = fits.open("/data/amarchal/ROHSA_paper/data/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2.fits") 
hdr = hdu[0].header

cube = hdu[0].data
cube[np.where(cube != cube)] = 0.

dim_x = cube.shape[2]
dim_y = cube.shape[1]
dim_v = cube.shape[0]

CDELT = hdr['CDELT3'] #*1.e-3  #km.s-1
CRVAL = hdr['CRVAL3'] #*1.e-3  #km.s-1
CRPIX = hdr['CRPIX3']

reso = np.abs(CDELT)

#ROHSA descomposition
N_G = 8
cube_all = []
for idx in np.arange(N_G):
    data = np.genfromtxt("/data/amarchal/ROHSA_paper/ROHSA/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_{}.dat".format(idx)) 
    
    id1 = data[:, 0]
    id2 = data[:, 1]
    sigma = data[:, 4] * reso
    amp = data[:, 2]
    mean = data[:, 3] - 1
    
    reconstructed_cube = np.zeros((dim_v, dim_y, dim_x))
    
    for i in range(len(amp)):
        gauss = gaussian(np.arange(dim_v), amp[i], mean[i], sigma[i]/reso)
        reconstructed_cube[:, int(id1[i]), int(id2[i])] += gauss
    
    n_gauss = len(id1) / (dim_y*dim_x)
    params = np.zeros((3*n_gauss, dim_y, dim_x))

    v = mean2vel(CRVAL, CDELT, CRPIX, np.arange(cube.shape[0]))
    vmean = mean2vel(CRVAL, CDELT, CRPIX, mean)

    i__ = 0
    for i in range(dim_y):
        for j in range(dim_x):
            for k in range(n_gauss):
                params[0+(3*k),i,j] = amp[i__]
                params[1+(3*k),i,j] = vmean[i__]
                params[2+(3*k),i,j] = sigma[i__]
                i__ += 1

    fields = [np.sqrt(2.*np.pi) * params[0+(3*k)] * params[2+(3*k)] for k in np.arange(n_gauss)]
    field = [f * C / 1.e18 for f in fields]
    ampfield = [params[0+(3*k)] for k in np.arange(n_gauss)]
    vfield = [params[1+(3*k)] for k in np.arange(n_gauss)]
    sigfield = [params[2+(3*k)] for k in np.arange(n_gauss)]
    
    cube_all.append(reconstructed_cube)

residual = [(cube - rohsa_cube) / cube for rohsa_cube in cube_all]
residual_NHI = [(np.sum(cube,0) - np.sum(rohsa_cube,0)) / np.sum(cube,0) for rohsa_cube in cube_all]
mu3 = [np.abs(moment((cube - foo).ravel(), moment=3, nan_policy='propagate')) for foo in cube_all]
mu3_NHI = [np.abs(moment((np.sum(cube,0) - np.sum(foo,0)).ravel(), moment=3, nan_policy='propagate')) for foo in cube_all]

fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.tick_params(labelsize=14)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
bins = np.linspace(-10, 20, 100)
ax.set_xlim([-10., 20])
ax.set_xlabel(r'(N$_{HI}$ - $\tilde{N}_{HI}$) / N$_{HI}$ [%]' ,  fontsize = 16)
ax.set_ylabel(r'Normalized PDF',  fontsize = 16)
ax.hist(residual_NHI[1].ravel()*100., bins=bins, histtype="step", normed=True, linewidth=2.5, label=r"N = 1; $|\mu_3|$ = {:{width}.{prec}f}".format(mu3_NHI[1],width=5,prec=2))
ax.hist(residual_NHI[3].ravel()*100., bins=bins, histtype="step", normed=True, linewidth=2.5, label=r"N = 3; $|\mu_3|$ = {:{width}.{prec}f}".format(mu3_NHI[3],width=5,prec=2))
ax.hist(residual_NHI[5].ravel()*100., bins=bins, histtype="step", normed=True, linewidth=2.5, label=r"N = 5; $|\mu_3|$ = {:{width}.{prec}f}".format(mu3_NHI[5],width=5,prec=2))
ax.hist(residual_NHI[7].ravel()*100., bins=bins, histtype="step", normed=True, linewidth=2.5, label=r"N = 7; $|\mu_3|$ = {:{width}.{prec}f}".format(mu3_NHI[7],width=5,prec=2))
ax.hist(residual_NHI[0].ravel()*100., bins=bins, histtype="step", normed=True, linewidth=2.5, label=r"N = 8; $|\mu_3|$ = {:{width}.{prec}f}".format(mu3_NHI[0],width=5,prec=2))
ax.legend()
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
plt.setp(ax.get_yticklabels()[0], visible=False)   
plt.setp(ax.get_xticklabels()[0], visible=False)   
plt.setp(ax.get_xticklabels()[-1], visible=False)   
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.legend(loc = 1, numpoints = 1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize = 'small')
plt.savefig('plot/PDF_resiudal_SIMU.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)
