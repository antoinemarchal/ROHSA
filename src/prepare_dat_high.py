import numpy as np
from astropy.io import fits
from astropy import constants as const
from astropy import units as u

from ROHSApy import ROHSA
import marchalib as ml

#Open data GNILC_IRIS-SFD_adim
fitsname = "/mnt/raid-cita/amarchal/DUST/data/CIBA/intensity_adim.fits"
hdu = fits.open(fitsname)
hdr = hdu[0].header
cube = hdu[0].data[:,:256,:256] *1.e6

fitsname_noise = "/mnt/raid-cita/amarchal/DUST/data/CIBA/intensity_adim_noise.fits"
hdu_noise = fits.open(fitsname_noise)
hdr_noise = hdu_noise[0].header
rms_cube = hdu_noise[0].data[:,:256,:256]  *1.e6
# rms_cube = np.full((cube.shape[0],cube.shape[1],cube.shape[2]), 1.)

fitsname_hi = "/mnt/raid-cita/amarchal/DUST/data/CIBA/NHI_HI4PI_glon140_glat85.fits"
hdu_hi = fits.open(fitsname_hi)
hdr_hi = hdu_hi[0].header
cube_hi = hdu_hi[0].data[:256,:256]
NHI = np.ones((2,cube.shape[1],cube.shape[2]))
NHI[1] = cube_hi

#Open color correction file
color = fits.open("/home/amarchal/ROHSA/data/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

core_cube = ROHSA(cube)
core_NHI = ROHSA(NHI)
core_color = ROHSA(np.zeros((0,color.shape[0],color.shape[1])))
core_rms = ROHSA(rms_cube)

core_cube.cube2dat(filename="/home/amarchal/ROHSA/src/SED.dat")
core_NHI.cube2dat(filename="/home/amarchal/ROHSA/src/NHI.dat")
core_rms.cube2dat(filename="/home/amarchal/ROHSA/src/rms.dat")
core_color.rms_map(color, filename="/home/amarchal/ROHSA/src/color.dat")
