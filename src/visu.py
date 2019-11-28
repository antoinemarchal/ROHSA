import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy import constants as const
from astropy import units as u
from scipy.special import gamma

from ROHSApy import ROHSA
import turbulence as tb

import mod_optimize as mod_opt
reload(mod_opt)

def planckct():
    colombi1_cmap = matplotlib.colors.ListedColormap(np.loadtxt("/home/amarchal/library/magnetar/Planck_Parchment_RGB.txt")/255.)
    colombi1_cmap.set_bad("white")
    colombi1_cmap.set_under("white")
    return colombi1_cmap


def poly(coef, xx, yy, degree):
    model = 0.
    k = 0
    for i in np.arange(degree):
        for j in np.arange(degree+1-i):
            model += coef[k] * xx**i * yy**j
            k += 1
    return model

#Open data GNILC_IRIS-SFD_adim
path = "/data/amarchal/GNILC_IRIS-SFD_adim/"
# fitsname = "GNILC_IRIS-SFD_adim_G86.fits"
fitsname = "GNILC_IRIS-SFD_adim_N1.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
# cube = hdu[0].data[:,44-32:44+32,44-32:44+32] *1.e6   #Attention rescale
cube = hdu[0].data[:,32-32:32+32,54-32:54+32] *1.e6   #Attention rescale N1

freq = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.array([353,545,857,3000])])*u.GHz
wavelength = (const.c / freq).to(u.micron)
wavelength_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

#Open color correction file
color = fits.open("/data/amarchal/PLANCK/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

#Open HI model CNM WNM with ROHSA
path_HI = "/data/amarchal/EN/model/"
fitsname_HI = "GHIGLS_N1_LVC_HVC.fits"
hdu_HI = fits.open(path_HI+fitsname_HI)
hdr_HI = hdu_HI[0].header
NHI = hdu_HI[0].data  *u.cm**-2 # 1.e18 

# NHI = np.array([NHI[0].value])  *u.cm**-2

# NHI = np.zeros((1,cube.shape[1],cube.shape[2]))
# NHI[0,:,:] = np.sum(hdu_HI[0].data[:,44-32:44+32,44-32:44+32],0)
# NHI = NHI*u.cm**-2

rms_map = np.ones((cube.shape[1],cube.shape[2]))

core_cube = ROHSA(cube)
core_NHI = ROHSA(NHI.value)

gaussian = core_cube.read_gaussian("./SED_mbb_run_0.dat")
gaussian[1::3] += 1. #ATTENtioN beta -1 ICI car function gauss used with ROHSApy

wavelength_mean_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

model_full = 0.
model = 0.
l0 = 349.81617036*u.micron
for k in np.arange(NHI.shape[0]):
    cc = np.array([poly(color[i,:],gaussian[1+(k*3)],gaussian[2+(k*3)],5) for i in np.arange(len(freq[:,0,0]))])
    model_full += mod_opt.MBB_l_adim(wavelength_mean_full,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0)
    model += mod_opt.MBB_l_adim(wavelength,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0) 

mbbcc = model / cc #FIXME model * cc

sigfield = gaussian[0::3]
betafield = gaussian[1::3]
Tfield = gaussian[2::3]

luminosity = gaussian[0::3]*(gaussian[2::3]**(4.+gaussian[1::3]) * gamma(4.+gaussian[1::3]) * (const.k_B.cgs*l0.cgs/const.h.cgs/const.c.cgs).value**gaussian[1::3])

stop

#PLOT SIGMA
lh = 2; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10,12.1))
fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im = axs[i][j].imshow(sigfield[k], origin="lower", cmap=planckct())
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend="both")
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=14.) 
        if i == lh-1 : cbar.set_label(r"$\sigma$ / [arb. unit]", fontsize=16.)
        k += 1
plt.savefig('plot/mosaic_sigfield.pdf', format='pdf')

#PLOT BETA
lh = 2; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10,12.1))
fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im = axs[i][j].imshow(betafield[k], origin="lower", cmap=planckct())
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend="both")
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=14.) 
        if i == lh-1 : cbar.set_label(r"$\beta$ / [arb. unit]", fontsize=16.)
        k += 1
plt.savefig('plot/mosaic_betafield.pdf', format='pdf')

#PLOT T
lh = 2; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10,12.1))
fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im = axs[i][j].imshow(Tfield[k], origin="lower", cmap=planckct())
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend="both")
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=14.) 
        if i == lh-1 : cbar.set_label(r"T / [arb. unit]", fontsize=16.)
        k += 1
plt.savefig('plot/mosaic_Tfield.pdf', format='pdf')
