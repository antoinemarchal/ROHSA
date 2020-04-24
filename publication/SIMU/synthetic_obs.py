#!/home/amarchal/py2env/bin/python
'''This program build synthetic obs (21cm line) from T,n and vz which are the three-dimensional 
field of the numerical simulation based on the work of Saury et al. 2014'''

import numpy as np
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units
from astropy import constants as const
from scipy import ndimage
import scipy.integrate as integrate
import FITS_tools

plt.ion()

plot = False
cm = plt.get_cmap('viridis')
cm.set_bad(color='black')
imkw = dict(origin='lower', interpolation='none', cmap=cm)

def I_Tb(params):
    (u, vz, n_Delta, T, C, Delta2, dz) = params
    dI = n_Delta * np.exp(- (u - (vz))**2 / (2.*Delta2))
    dI[np.where(dI != dI)] = 0.
    I = 1./(C * np.sqrt(2.*np.pi)) * integrate.simps(dI, dx=dz, axis=0)
    return I

# Constant
m_h = 1.6737236e-27 #kg
C = 1.82243e18 #K-1cm-2 / (km.s-1)
pc2cm = units.pc.to(units.m) * 1.e2
box_size = 40. # pc
resolution = 1024.
dz = (box_size / resolution) * pc2cm

# Open data 
path_simu = '/data/amarchal/ROHSA_paper/data/Saury2014/'
path_out = '/data/amarchal/ROHSA_paper/data/synthetic_obs/'

hdu_list_rho = fits.open(path_simu + 'rho_016_subgrid_256.fits')
hdu_list_T = fits.open(path_simu + 'T_016_subgrid_256.fits')
hdu_list_vz = fits.open(path_simu + 'vz_016_subgrid_256.fits')

reso = 0.8 #km.s-1

rho_cube = hdu_list_rho[0].data #g.cm-3
T_cube = hdu_list_T[0].data
vz_cube = hdu_list_vz[0].data * 1.e-5 #km.s-1 ATTENTION 

## CUT TEMPERATURE
Tk_lim_inf = 0
Tk_lim_sup = np.inf

idx_phase = np.where((T_cube > Tk_lim_inf) & (T_cube < Tk_lim_sup))

rho_cube_phase = np.zeros((rho_cube.shape[0], rho_cube.shape[1], rho_cube.shape[2]))
T_cube_phase = np.zeros((rho_cube.shape[0], rho_cube.shape[1], rho_cube.shape[2]))
vz_cube_phase = np.zeros((rho_cube.shape[0], rho_cube.shape[1], rho_cube.shape[2]))

rho_cube_phase[idx_phase] = rho_cube[idx_phase]
T_cube_phase[idx_phase] = T_cube[idx_phase]
vz_cube_phase[idx_phase] = vz_cube[idx_phase]
##

# Preliminary calculation
Delta2 = ((const.k_B.value * T_cube_phase / m_h)) * 1.e-6 #km.s-1
n = rho_cube_phase/(m_h*1.e3)
n_Delta = n / np.sqrt(Delta2)

# Spectral range
u = np.arange(-40,40+reso, reso)

map_u = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))
for i in np.arange(T_cube_phase.shape[1]):
    for j in np.arange(T_cube_phase.shape[2]):
        map_u[:,i,j] = u

Tb = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))
Tb_thin = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))
tau_in_front = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))

for i in tqdm(range(T_cube_phase.shape[0])):    
    Tb_z = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))
    tau_z = 1. / (C * np.sqrt(2.*np.pi)) * n_Delta[i] / T_cube_phase[i] * np.exp(- (map_u - (vz_cube_phase[i]))**2 / (2.*Delta2[i])) * dz
    idx_nonzero = ~np.isnan(tau_z[0])
    Tb_z[:,idx_nonzero] = T_cube_phase[i,idx_nonzero] * (1. - np.exp(-1.*tau_z[:,idx_nonzero])) * np.exp(-1.*tau_in_front[:,idx_nonzero])
    tau_in_front[:,idx_nonzero] += tau_z[:,idx_nonzero]
    Tb += Tb_z
    Tb_thin[:,idx_nonzero] += tau_z[:,idx_nonzero] * T_cube_phase[i,idx_nonzero]
 
# Tb_thin_fast = np.zeros((len(u), T_cube_phase.shape[1], T_cube_phase.shape[2]))
# for i in tqdm(range(len(u))):
#     Tb_thin_fast[i] = I_Tb((u[i], vz_cube_phase, n_Delta, T_cube_phase, C, Delta2, dz))    

fileout = 'Tb_reso_' + str(reso) + 'km.s-1_' + "Tmin_" + str(Tk_lim_inf) + "_Tmax_" + str(Tk_lim_sup) + '_ROHSA.fits'
fileout_thin = 'Tb_reso_' + str(reso) + 'km.s-1_' + "Tmin_" + str(Tk_lim_inf) + "_Tmax_" + str(Tk_lim_sup) + '_ROHSA_thin.fits'

# Write PPV cube
hdu0 = fits.PrimaryHDU(Tb)
hdu0.header['COMMENT'] = 'Brightness Temperature Tb'
hdu0.header['NAXIS'] = 3
hdu0.header['NAXIS1'] = Tb.shape[1]
hdu0.header['NAXIS2'] = Tb.shape[2]
hdu0.header['NAXIS3'] = len(u)
hdu0.header['CTYPE3'] = 'v [km.s-1]'
hdu0.header['CRVAL3'] = u[40]
hdu0.header['CDELT3'] = reso
hdu0.header['CRPIX3'] = 40
hdu0.header['BUNIT'] = 'K'
hdulist = fits.HDUList([hdu0])
hdulist.writeto(path_out + fileout, overwrite=True)
    
# Write PPV cube thin limit
hdu0 = fits.PrimaryHDU(Tb_thin)
hdu0.header['COMMENT'] = 'Brightness Temperature Tb'
hdu0.header['NAXIS'] = 3
hdu0.header['NAXIS1'] = Tb_thin.shape[1]
hdu0.header['NAXIS2'] = Tb_thin.shape[2]
hdu0.header['NAXIS3'] = len(u)
hdu0.header['CTYPE3'] = 'v [km.s-1]'
hdu0.header['CRVAL3'] = u[40]
hdu0.header['CDELT3'] = reso
hdu0.header['CRPIX3'] = 40
hdu0.header['BUNIT'] = 'K'
hdulist = fits.HDUList([hdu0])
hdulist.writeto(path_out + fileout_thin, overwrite=True)
