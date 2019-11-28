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
path = "/data/amarchal/TAURUS/data/"
hdu = fits.open(path+"Planck_PR3_glon170_glat-16_adim.fits")
cube = hdu[0].data *1.e4
hdr = hdu[0].header
w = tb.proj.wcs2D(hdr)

cube_noise = fits.open(path+"Planck_PR3_glon170_glat-16_adim_noise.fits")[0].data *1.e4

core_cube = ROHSA(cube)

NHI = np.ones((1,cube.shape[1],cube.shape[2]))

freq = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.array([353,545,857,3000])])*u.GHz
wavelength = (const.c / freq).to(u.micron)
wavelength_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

#Open color correction file
color = fits.open("/data/amarchal/PLANCK/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

gaussian = core_cube.read_gaussian("./SED_mbb_run_stefan.dat")
gaussian[1::3] += 1. #ATTENtioN beta -1 ICI car function gauss used with ROHSApy

wavelength_mean_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

model_full = 0.
model = 0.
l0 = 349.81617036*u.micron
for k in np.arange(NHI.shape[0]):
    cc = np.array([poly(color[i,:],gaussian[1+(k*3)],gaussian[2+(k*3)],5) for i in np.arange(len(freq[:,0,0]))])
    model_full += mod_opt.MBB_l_adim(wavelength_mean_full,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0)
    model += mod_opt.MBB_l_adim(wavelength,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0) 

mbbcc = model / cc

sigfield = gaussian[0::3]
betafield = gaussian[1::3]
Tfield = gaussian[2::3]

luminosity = gaussian[0::3]*(gaussian[2::3]**(4.+gaussian[1::3]) * gamma(4.+gaussian[1::3]) * (const.k_B.cgs*l0.cgs/const.h.cgs/const.c.cgs).value**gaussian[1::3])

#Plot integrated column density field TOT
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(np.log10(sigfield[0]), origin="lower", cmap="viridis")
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"$log_{10}(\tau_{350 \mu m})$ / [arb. units]", fontsize=18.)
plt.savefig("plot/" + 'sig_taurus.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot integrated column density field TOT
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(betafield[0], vmin=1.3, **tb.imkw_cubehelix)
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"$\beta$", fontsize=18.)
plt.savefig("plot/" + 'beta_taurus.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot integrated column density field TOT
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(Tfield[0], origin="lower", cmap="cividis")
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"$T_{dust}$", fontsize=18.)
plt.savefig("plot/" + 'Tdust_taurus.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)


