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


#Open data GNILC_IRIS-SFD_adim
fitsname = "/data/amarchal/SIMUSED/data/intensity_adim.fits"
hdu = fits.open(fitsname)
hdr = hdu[0].header
cube = hdu[0].data

core_cube = ROHSA(cube)

NHI = np.ones((1,cube.shape[1],cube.shape[2]))

freq = np.array([np.full((cube.shape[0],cube.shape[1]),i) for i in np.array([3000,857,545,353])])*u.GHz
wavelength = (const.c / freq).to(u.micron)
wavelength_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

#Open color correction file
color = fits.open("/data/amarchal/PLANCK/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

gaussian = core_cube.read_gaussian("./SED_mbb_run_0.dat")
gaussian[1::3] += 1. #ATTENTION beta -1

stop

wavelength_mean_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

model_full = 0.
model = 0.
l0 = 349.81617036*u.micron
for k in np.arange(NHI.shape[0]):
    cc = np.array([poly(color[i,:],gaussian[1+(k*3)],gaussian[2+(k*3)],5) for i in np.arange(len(freq[:,0,0]))])
    model_full += mod_opt.MBB_l_adim(wavelength_mean_full,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0)
    model += mod_opt.MBB_l_adim(wavelength,NHI[k],gaussian[0+(k*3)]*u.cm**2,gaussian[1+(k*3)],gaussian[2+(k*3)]*u.K,l0) 

sigfield = gaussian[0::3]
betafield = gaussian[1::3]
Tfield = gaussian[2::3]
