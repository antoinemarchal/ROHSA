import numpy as np
from astropy.io import fits
from astropy import constants as const
from astropy import units as u

from ROHSApy import ROHSA
import turbulence as tb

#Open data GNILC_IRIS-SFD_adim
path = "/data/amarchal/TAURUS/data/"
cube = fits.open(path+"Planck_PR3_glon170_glat-16_adim.fits")[0].data *1.e4
cube_noise = fits.open(path+"Planck_PR3_glon170_glat-16_adim_noise.fits")[0].data *1.e4

NHI = np.ones((1,cube.shape[1],cube.shape[2]))

# rms_cube = np.ones((cube.shape[0],cube.shape[1],cube.shape[2]))

#Open color correction file
color = fits.open("/data/amarchal/PLANCK/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

core_cube = ROHSA(cube)
core_NHI = ROHSA(NHI)
core_color = ROHSA(np.zeros((0,color.shape[0],color.shape[1])))
core_rms = ROHSA(cube_noise)
# core_rms = ROHSA(rms_cube)

core_cube.cube2dat(filename="/home/amarchal/ROHSA/src/SED.dat")
core_NHI.cube2dat(filename="/home/amarchal/ROHSA/src/NHI.dat")
core_rms.cube2dat(filename="/home/amarchal/ROHSA/src/rms.dat")
core_color.rms_map(color, filename="/home/amarchal/ROHSA/src/color.dat")

