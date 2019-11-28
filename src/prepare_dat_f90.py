import numpy as np
from astropy.io import fits
from astropy import constants as const
from astropy import units as u

from ROHSApy import ROHSA
import turbulence as tb

#Open data GNILC_IRIS-SFD_adim
path = "/data/amarchal/GNILC_IRIS-SFD_adim/"
fitsname = "GNILC_IRIS-SFD_adim_N1.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
# cube = hdu[0].data[:,44-32:44+32,44-32:44+32] *1.e6   #Attention rescale G86
cube = hdu[0].data[:,32-32:32+32,54-32:54+32] *1.e6   #Attention rescale N1

freq = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.array([353,545,857,3000])])*u.GHz
wavelength = (const.c / freq).to(u.micron)
wavelength_full = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.arange(100,900,1)]) * u.micron

#Open HI model CNM WNM with ROHSA
path_HI = "/data/amarchal/EN/model/"
fitsname_HI = "GHIGLS_N1_LVC_HVC.fits"
hdu_HI = fits.open(path_HI+fitsname_HI)
hdr_HI = hdu_HI[0].header
NHI = hdu_HI[0].data *u.cm**-2 # 1.e18 

# NHI = np.array([NHI[0].value])  *u.cm**-2

# NHI = np.zeros((1,cube.shape[1],cube.shape[2]))
# NHI[0,:,:] = np.sum(hdu_HI[0].data[:,44-32:44+32,44-32:44+32],0)
# NHI = NHI*u.cm**-2

rms_cube = np.ones((cube.shape[0],cube.shape[1],cube.shape[2]))

#Open color correction file
color = fits.open("/data/amarchal/PLANCK/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

core_cube = ROHSA(cube)
core_NHI = ROHSA(NHI.value)
core_color = ROHSA(np.zeros((0,color.shape[0],color.shape[1])))
core_rms = ROHSA(rms_cube)

core_cube.cube2dat(filename="/home/amarchal/ROHSA/src/SED.dat")
core_NHI.cube2dat(filename="/home/amarchal/ROHSA/src/NHI.dat")
core_rms.cube2dat(filename="/home/amarchal/ROHSA/src/rms.dat")
core_color.rms_map(color, filename="/home/amarchal/ROHSA/src/color.dat")


