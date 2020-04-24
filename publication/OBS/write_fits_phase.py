import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import turbulence as tb
from ROHSApy import ROHSA

#DATA                                                                                                                                                                                                                                                                                 
hdu = fits.open("/data/amarchal/ROHSA_paper/data/observation/GHIGLS_NEP_Tb.fits")
hdr = hdu[0].header
hdr['CRPIX3'] -= 200; hdr['NAXIS3'] = 300
w = tb.proj.wcs2D(hdr)
cube = hdu[0].data[0][200:500,:,:]
core = ROHSA(cube, hdr=hdr)

#FULL FIELD                                                                                                                                                                                                                                                                           
core = ROHSA(np.zeros((cube.shape[0],cube.shape[1],cube.shape[2])), hdr=hdr)

output = core.read_gaussian("/data/amarchal/ROHSA_paper/ROHSA/GHIGLS_NEP_Tb_gauss_run_4.dat")

params = core.physical_gaussian(output)
reconstructed_cube = core.return_result_cube(output)

field = params[0::3] * params[2::3] * np.sqrt(2.*np.pi) * tb.cst.C / 1.e18
ampfield = params[0::3]
vfield = params[1::3]
sigfield = params[2::3]

ampfield_pix = output[0::3]
vfield_pix = output[1::3]
sigfield_pix = output[2::3]

iddx = np.argsort(np.mean(vfield, axis=(1,2)))

field = [field[idd] for idd in iddx]
vfield = [vfield[idd] for idd in iddx]
ampfield = [ampfield[idd] for idd in iddx]
sigfield = [sigfield[idd] for idd in iddx]

vfield_pix = [vfield_pix[idd] for idd in iddx]
ampfield_pix = [ampfield_pix[idd] for idd in iddx]
sigfield_pix = [sigfield_pix[idd] for idd in iddx]

#CNM/LNM/WNM fraction                                                                                                                                                                                                                                                                 
#Sort sigfield                                                                                                                                                                                                                                                                        
iddx = np.argsort(np.median(sigfield, axis=(1,2)))

field = [field[idd] for idd in iddx]
vfield = [vfield[idd] for idd in iddx]
ampfield = [ampfield[idd] for idd in iddx]
sigfield = [sigfield[idd] for idd in iddx]

vfield_pix = [vfield_pix[idd] for idd in iddx]
ampfield_pix = [ampfield_pix[idd] for idd in iddx]
sigfield_pix = [sigfield_pix[idd] for idd in iddx]

idx_CNM_local = np.where((np.median(sigfield, axis=(1,2)) < 3.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_LNM_local = np.where((np.median(sigfield, axis=(1,2)) > 3.) & (np.median(sigfield, axis=(1,2)) < 6.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_WNM_local = np.where((np.median(sigfield, axis=(1,2)) > 6.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]

idx_CNM_ivc = np.where((np.median(sigfield, axis=(1,2)) < 3.) & (np.median(vfield, axis=(1,2)) < -20.) | (np.median(vfield, axis=(1,2)) > 20.))[0]
idx_LNM_ivc = np.where((np.median(sigfield, axis=(1,2)) > 3.) & (np.median(sigfield, axis=(1,2)) < 8.) & (np.median(vfield, axis=(1,2)) < -20.) | (np.median(vfield, axis=(1,2)) > 20.))[0]
idx_WNM_ivc = np.where((np.median(sigfield, axis=(1,2)) > 8.) & (np.median(vfield, axis=(1,2)) < -20.) | (np.median(vfield, axis=(1,2)) > 20.))[0]

#Reconstruct cube phase                                                                                                                                                                                                                                                               
model_CNM_local = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_CNM_local], pixfield=np.array(vfield_pix)[idx_CNM_local], sigfield=np.array(sigfield_pix)[idx_CNM_local])
model_WNM_local = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_WNM_local], pixfield=np.array(vfield_pix)[idx_WNM_local], sigfield=np.array(sigfield_pix)[idx_WNM_local])
model_LNM_local = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_LNM_local], pixfield=np.array(vfield_pix)[idx_LNM_local], sigfield=np.array(sigfield_pix)[idx_LNM_local])

model_CNM_ivc = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_CNM_ivc], pixfield=np.array(vfield_pix)[idx_CNM_ivc], sigfield=np.array(sigfield_pix)[idx_CNM_ivc])
model_WNM_ivc = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_WNM_ivc], pixfield=np.array(vfield_pix)[idx_WNM_ivc], sigfield=np.array(sigfield_pix)[idx_WNM_ivc])
model_LNM_ivc = core.return_result_cube(ampfield=np.array(ampfield_pix)[idx_LNM_ivc], pixfield=np.array(vfield_pix)[idx_LNM_ivc], sigfield=np.array(sigfield_pix)[idx_LNM_ivc])

stop

#Write output local                                                                                                                                                                                                                                                                        
hdu0 = fits.PrimaryHDU(model_CNM_local)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_CNM_local.fits", overwrite=True)

hdu0 = fits.PrimaryHDU(model_LNM_local)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_LNM_local.fits", overwrite=True)

hdu0 = fits.PrimaryHDU(model_WNM_local)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_WNM_local.fits", overwrite=True)

#Write output ivc                                                                                                                                                                                                                                                                        
hdu0 = fits.PrimaryHDU(model_CNM_ivc)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_CNM_ivc.fits", overwrite=True)

hdu0 = fits.PrimaryHDU(model_LNM_ivc)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_LNM_ivc.fits", overwrite=True)

hdu0 = fits.PrimaryHDU(model_WNM_ivc)
hdu0.header = hdr
hdulist = fits.HDUList([hdu0])
hdulist.writeto("/home/amarchal/Projects/ROHSA_paper/OBS/output/GHIGLS_NEP_Tb_WNM_ivc.fits", overwrite=True)
