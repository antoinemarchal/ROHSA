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
import colorcet as cc
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
import matplotlib.axes as maxes

import turbulence as tb

plt.ion()

# cm_coolwarm = plt.get_cmap('pink')
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
hdu = fits.open("/data/amarchal/ROHSA_paper/data/observation/GHIGLS_NEP_Tb.fits") 
data = np.genfromtxt("/data/amarchal/ROHSA_paper/ROHSA/GHIGLS_NEP_Tb_gauss_run_4.dat") # 4

cube = hdu[0].data[0]
cube[np.where(cube != cube)] = 0.

cube = cube[200:500,:,:]

hdr = hdu[0].header

dim_x = cube.shape[2]
dim_y = cube.shape[1]
dim_v = cube.shape[0]

id1 = data[:, 0]
id2 = data[:, 1]
sigma = data[:, 4]
amp = data[:, 2]
mean = data[:, 3] - 1

CDELT = hdr['CDELT3'] *1.e-3  #km.s-1  
CRVAL = hdr['CRVAL3'] *1.e-3  #km.s-1                                                                                   
CRPIX = hdr['CRPIX3'] - 200 

reso = np.abs(CDELT)

v = mean2vel(CRVAL, CDELT, CRPIX, np.arange(cube.shape[0]))
vmean = mean2vel(CRVAL, CDELT, CRPIX, mean)

rms = np.zeros((cube.shape[1],cube.shape[2]))
for i in np.arange(cube.shape[1]):
    for j in np.arange(cube.shape[2]):
        rms[i,j] = np.std(cube[:40,i,j])

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

original_map = np.sum(cube, axis=0) * reso * C / 1.e18 / 100. #10^20
reconstructed_map = np.sum(reconstructed_cube, axis=0) * C * reso / 1.e18 / 100. #10^20

mu3 = np.abs(moment((reconstructed_cube - cube).ravel(), moment=3, nan_policy='propagate'))

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

mean_sig = [np.around(np.mean(sigfields), decimals=1) for sigfields in sigfield]
mean_v = [np.around(np.mean(vfields), decimals=1) for vfields in vfield]

stop

#HEATMAP
x_bins = np.linspace(np.min(vmean), np.max(vmean), 1200)
# y_bins = np.logspace(np.log(0.1), np.log(np.max(sigma)), 1200)
y_bins = np.linspace(0.1, np.max(sigma), 1200)
H, xedges, yedges = np.histogram2d(vmean, sigma*reso, weights=np.sqrt(2.*np.pi)*amp*(sigma*reso)*C/1.e18/np.sum(field), bins=[x_bins, y_bins])
H = np.ma.masked_invalid(np.atleast_2d(H))
fig = plt.figure(figsize=(16.,8.))
ax = fig.add_subplot(111)
# ax.set_yscale('log')
# ax.set_ylim([0.8, 30.])
ax.set_xlim([np.min(vmean),np.max(vmean)])
ax.set_xlabel(r'v [km s$^{-1}$]',  fontsize = 16)
ax.set_ylabel(r'$\sigma$ [km s$^{-1}$]',  fontsize = 16)
ax.yaxis.set_major_formatter(ScalarFormatter()) 
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.05)
im = ax.pcolormesh(xedges, yedges, np.log10(H.T), cmap=cm_inferno, vmin=-8., vmax=-2.5)
cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.set_label(r'log$_{10}$($N_{HI}$) [fraction of total emission]', fontsize = 16)
plt.savefig("plot/" + 'heatmap_NEP_linear.png', format='png', bbox_inches='tight', pad_inches=0.02)

#CNM/LNM/WNM fraction
#Sort sigfield
iddx = np.argsort(np.median(sigfield, axis=(1,2)))

field = [field[idd] for idd in iddx]
vfield = [vfield[idd] for idd in iddx]
ampfield = [ampfield[idd] for idd in iddx]
sigfield = [sigfield[idd] for idd in iddx]

idx_CNM_local = np.where((np.median(sigfield, axis=(1,2)) < 3.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_LNM_local = np.where((np.median(sigfield, axis=(1,2)) > 3.) & (np.median(sigfield, axis=(1,2)) < 6.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_WNM_local = np.where((np.median(sigfield, axis=(1,2)) > 6.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_INTER_local = np.where((np.median(sigfield, axis=(1,2)) > 3.) & (np.median(vfield, axis=(1,2)) > -20.) & (np.median(vfield, axis=(1,2)) < 20.))[0]
idx_CNM_ivc = np.where((np.median(sigfield, axis=(1,2)) < 3.) & (np.median(vfield, axis=(1,2)) < -20.))[0]
idx_LNM_ivc = np.where((np.median(sigfield, axis=(1,2)) > 3.) & (np.median(sigfield, axis=(1,2)) < 8.) & (np.median(vfield, axis=(1,2)) < -20.))[0]
idx_WNM_ivc = np.where((np.median(sigfield, axis=(1,2)) > 8.) & (np.median(vfield, axis=(1,2)) < -20.))[0]

NHI_CNM_local = sum(field[i] for i in idx_CNM_local)
NHI_LNM_local = sum(field[i] for i in idx_LNM_local)
NHI_WNM_local = sum(field[i] for i in idx_WNM_local)
NHI_INTER_local = NHI_LNM_local + NHI_WNM_local
NHI_TOT_local = NHI_CNM_local + NHI_INTER_local
F_CNM_local = NHI_CNM_local / (NHI_INTER_local)

NHI_CNM_ivc = sum(field[i] for i in idx_CNM_ivc)
NHI_LNM_ivc = sum(field[i] for i in idx_LNM_ivc)
NHI_WNM_ivc = sum(field[i] for i in idx_WNM_ivc)

NHI_TOT_ivc = NHI_CNM_ivc+NHI_LNM_ivc+NHI_WNM_ivc

cube_INTER_ROHSA = np.zeros(cube.shape)
for k in idx_INTER_local:
    for i in np.arange(cube.shape[1]):
        for j in np.arange(cube.shape[2]):
            cube_INTER_ROHSA[:,i,j] += gaussian(np.arange(dim_v), ampfield[k][i,j], vel2mean(CRVAL, CDELT, CRPIX, [vfield[k][i,j]])[0], sigfield[k][i,j]/reso) 

cube_CNM_ROHSA = np.zeros(cube.shape)
for k in idx_CNM_local:
    for i in np.arange(cube.shape[1]):
        for j in np.arange(cube.shape[2]):
            cube_CNM_ROHSA[:,i,j] += gaussian(np.arange(dim_v), ampfield[k][i,j], vel2mean(CRVAL, CDELT, CRPIX, [vfield[k][i,j]])[0], sigfield[k][i,j]/reso)

CV_INTER_LOCAL_ROHSA = np.tensordot(v, cube_INTER_ROHSA, axes=([0],[0])) / np.sum(cube_INTER_ROHSA, axis=0)
DISP_INTER_LOCAL_ROHSA = np.sqrt(np.tensordot(np.array(v)**2, cube_INTER_ROHSA, axes=([0],[0])) / np.sum(cube_INTER_ROHSA, axis=0) - (CV_INTER_LOCAL_ROHSA)**2)

CV_CNM_LOCAL_ROHSA = np.tensordot(v, cube_CNM_ROHSA, axes=([0],[0])) / np.sum(cube_CNM_ROHSA, axis=0)
DISP_CNM_LOCAL_ROHSA = np.sqrt(np.tensordot(np.array(v)**2, cube_CNM_ROHSA, axes=([0],[0])) / np.sum(cube_CNM_ROHSA, axis=0) - (CV_CNM_LOCAL_ROHSA)**2)

#Iterpolation cv field where no data
#Convolution to interpolate missing data
CV_CNM_LOCAL_ROHSA[np.where(CV_CNM_LOCAL_ROHSA == 0.)] = np.nan
kernel = Gaussian2DKernel(0.5)
CV_CNM_LOCAL_ROHSA = convolve(CV_CNM_LOCAL_ROHSA, kernel)

stop

# #Open LIC NEP
# hdu_LIC = fits.open("NEP_LIC.fits") 
# LIC = hdu_LIC[0].data

# #POWER SPECTRUM NHI
# stat_cnm = tb.PowerS(NHI_LOCAL_CNM)
# stat_inter = tb.PowerS(NHI_LOCAL_INTER)
# stat_tot = tb.PowerS(NHI_LOCAL_CNM + NHI_LOCAL_INTER)
# sps1d_cnm, std_cnm = stat_cnm.sps1d(return_log=False)
# sps1d_inter, std_inter = stat_inter.sps1d(return_log=False)
# sps1d_tot, std_tot = stat_tot.sps1d(return_log=False)

# ks = stat_inter.get_ksps1d(unit_length=1)
# ksup = np.where(np.array(ks) > 0.4)[0][0]

# #POWER SPECTRUM NHI
# stat_cnm = tb.PowerS(vfield[9])
# stat_inter = tb.PowerS(CV_INTER_LOCAL_ROHSA)
# sps1d_cnm, std_cnm = stat_cnm.sps1d(return_log=False)
# sps1d_inter, std_inter = stat_inter.sps1d(return_log=False)

# ks = stat_inter.get_ksps1d(unit_length=1)
# ksup = np.where(np.array(ks) > 0.4)[0][0]

# PLOT HIST SIGMA  / AMP                                                                                                                
fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
bins = np.linspace(np.min(sigma*reso), np.max(sigma*reso), 800)
ax.hist(sigma*reso, weights=np.sqrt(2.*np.pi)*amp*(sigma*reso)/np.sum(field), bins=bins, log=False, histtype='step', color='k', normed=False, label=r'$\sigma weighted \, by \, A$', linewidth=2.5)
ax.set_xlim([0., 14.])
# ax.set_ylim([1., 3.75])  
ax.set_xlabel(r'$\sigma$ [km s$^{-1}$]', fontsize = 16)                                                                                                                                                                                            
ax.set_ylabel(r'PDF [fraction of total emission]', fontsize = 16)          
# ax.plot([9., 9.], [0., 4.], color="r", linestyle="--", linewidth=2.)
plt.savefig('plot/PDF_sigma_over_A_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02) 

# # PLOT HIST SIGMA  / AMP                                                                                                                
# idx_local = np.where((np.array(mean_v)> -15.) & (np.array(mean_v) < 15.))[0] 
# sig_local = np.array(list(sigfield[6].flatten()) + list(sigfield[7].flatten()) + list(sigfield[8].flatten()) + list(sigfield[9].flatten()) + list(sigfield[10].flatten()) + list(sigfield[11].flatten()))
# amp_local = np.array(list(field[6].flatten()) + list(field[7].flatten()) + list(field[8].flatten()) + list(field[9].flatten()) + list(field[10].flatten()) + list(field[11].flatten()))
# sig_nonlocal = np.array(list(sigfield[1].flatten()) + list(sigfield[2].flatten()) + list(sigfield[3].flatten()) + list(sigfield[4].flatten()) + list(sigfield[5].flatten()) + list(sigfield[12].flatten()))
# amp_nonlocal = np.array(list(field[1].flatten()) + list(field[2].flatten()) + list(field[3].flatten()) + list(field[4].flatten()) + list(field[5].flatten()) + list(field[12].flatten()))
# fig = plt.figure(figsize=(cm2inch((18.,18.))))
# ax = fig.add_subplot(111)
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.yaxis.set_major_formatter(ScalarFormatter())
# # bins = np.linspace(np.min(sigma*reso), np.max(sigma*reso), 800)
# bins_local = np.linspace(np.min(sig_local), np.max(sig_local), 800)
# bins_nonlocal = np.linspace(np.min(sig_nonlocal), np.max(sig_nonlocal), 800)
# ax.hist(sig_local, weights=amp_local/np.sum(field), bins=bins_local, log=False, histtype='step', color='r', normed=False, label=r'$\sigma weighted \, by \, A$', linewidth=2.5) #DEJA a sig reso dans field  donc OK !
# ax.hist(sig_nonlocal, weights=amp_nonlocal/np.sum(field), bins=bins_nonlocal, log=False, histtype='step', color='k', normed=False, label=r'$\sigma weighted \, by \, A$', linewidth=2.5) #DEJA a sig reso dans field  donc OK !
# ax.set_xlim([0., 17.])
# # ax.set_ylim([1., 3.75])  
# ax.set_xlabel(r'$\sigma$ [km s$^{-1}$]', fontsize = 16)                                                                                                                                                                                            
# ax.set_ylabel(r'PDF [fraction of total emission]', fontsize = 16)          
# # ax.plot([9., 9.], [0., 4.], color="r", linestyle="--", linewidth=2.)
# plt.savefig('plot/PDF_sigma_over_A_NEP_local.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02) 

w = wcs.WCS(naxis=2)
w.wcs.crpix = [hdr['CRPIX1'], hdr['CRPIX2']]
w.wcs.cdelt = np.array([hdr['CDELT1'], hdr['CDELT2']])
w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
w.wcs.ctype = [hdr['CTYPE1'], hdr['CTYPE2']]

#Plot integrated column density field TOT
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(original_map, aspect='auto', vmin=np.min(original_map), vmax=np.max(original_map), **imkw_inferno)
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)
plt.savefig("plot/" + 'NHI_TOT_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot integrated column density field TOT
fig = plt.figure(figsize=(18, 8))
# fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
ax0 = fig.add_subplot(131, projection=w)
ax0.set_xlabel(r"l", fontsize=16.)
ax0.set_ylabel(r"b", fontsize=16.)
im = ax0.imshow(original_map, vmin=np.min(original_map), vmax=np.max(original_map), **imkw_inferno)
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="2.5%", axes_class=maxes.Axes, pad=0.05)
cbar = fig.colorbar(im, cax=cax, format='%d', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=14.)
# cbar.set_label(r"N$_{HI}$ / [10$^{20}$ cm$^{-2}$]", fontsize=18.)

ax1 = fig.add_subplot(132, projection=w)
ax1.set_xlabel(r"l", fontsize=16.)
im = ax1.imshow(reconstructed_map, vmin=np.min(original_map), vmax=np.max(original_map), **imkw_inferno)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="2.5%", axes_class=maxes.Axes, pad=0.05)
cbar = fig.colorbar(im, cax=cax, format='%d', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=14.)

ax2 = fig.add_subplot(133, projection=w)
ax2.set_xlabel(r"l", fontsize=16.)
im = ax2.imshow(reconstructed_map-original_map, vmin=-0.1, vmax=0.1, **imkw_inferno)
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="2.5%", axes_class=maxes.Axes, pad=0.05)
cbar = fig.colorbar(im, cax=cax, format='%.2f', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=14.)
cbar.set_label(r"N$_{HI}$ / [10$^{20}$ cm$^{-2}$]", fontsize=16.)

plt.savefig("plot/" + 'NHI_TOT_NEP_model.png', format='png', bbox_inches='tight', pad_inches=0.02)

#Plot CNM FRACTION
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(F_CNM, aspect='auto', **imkw_inferno)
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"CNM fraction", fontsize=18.)
plt.savefig("plot/" + 'F_CNM_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

# #Plot integrated column density field and centroid velocity field LOCAL WNM
lh = 2; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10,12.))
fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)

im = axs[0][0].imshow(NHI_CNM_local, **imkw_inferno)
axs[0][0].axes.xaxis.set_ticklabels([])
axs[0][0].axis('off')
divider = make_axes_locatable(axs[0][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%d', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

im = axs[0][1].imshow(NHI_INTER_local, **imkw_inferno)
axs[0][1].axes.xaxis.set_ticklabels([])
axs[0][1].axis('off')
divider = make_axes_locatable(axs[0][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%d', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

im1 = axs[1][0].imshow(CV_CNM_LOCAL_ROHSA, **imkw_coolwarm)

axs[1][0].axes.xaxis.set_ticklabels([])
axs[1][0].axis('off')
divider = make_axes_locatable(axs[1][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=16.)

im1 = axs[1][1].imshow(CV_INTER_LOCAL_ROHSA, **imkw_coolwarm)
axs[1][1].axes.xaxis.set_ticklabels([])
axs[1][1].axis('off')
divider = make_axes_locatable(axs[1][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=16.)
plt.savefig('plot/NHI_CV_INTER_NEP.pdf', format='pdf')

#Plot mosaic field
lh = 4; lw = 3
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(cm2inch((22.,35.5))))
fig.subplots_adjust(top=1., bottom=0.05, left=0., right=1., hspace=0.12, wspace=0.01)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im = axs[i][j].imshow(field[k], **imkw_inferno)
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
        if i == lh-1 : cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=16.)
        k += 1
# plt.axis('off')
plt.savefig('plot/mosaic_field_NEP.pdf', format='pdf')

#Plot mosaic field mean
lh = 4; lw = 3
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(cm2inch((22.,35.5))))
fig.subplots_adjust(top=1., bottom=0.05, left=0., right=1., hspace=0.12, wspace=0.01)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        elev_min=np.nanmin(vfield[k])
        elev_max=np.nanmax(vfield[k])
        mid_val=0
        im1 = axs[i][j].imshow(vfield[k], **imkw_coolwarm)
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%d', extend="both")
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=14.) 
        if i == lh-1 : cbar.set_label(r"v [km s$^{-1}$]", fontsize=16.)
        k += 1
# plt.axis('off')
plt.savefig('plot/mosaic_vfield_NEP.pdf', format='pdf')

#Plot mosaic field sigma
lh = 4; lw = 3
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(cm2inch((22.,35.5))))
fig.subplots_adjust(top=1., bottom=0.05, left=0., right=1., hspace=0.12, wspace=0.01)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im1 = axs[i][j].imshow(sigfield[k], **imkw_cubehelix)
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend="both")
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=14.) 
        if i == lh-1 : cbar.set_label(r"$\sigma$ [km s$^{-1}$]", fontsize=16.)
        k += 1
# plt.axis('off')
plt.savefig('plot/mosaic_sigfield_NEP.pdf', format='pdf')

#Plot mosaic spectra                                                                                                                    
pvalues = np.logspace(-0.7, 0, 12)
pmin = pvalues[0]
pmax = pvalues[-1]

def norm(pval):
    return (pval - pmin) / float(pmax - pmin)

ny = 4; nx = 4                                                                                                                                          
center_y = 80; center_x = 80
cb = "magenta"
cw = "crimson"
fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(cm2inch((18.,14.))))
fig.subplots_adjust(hspace=0, wspace=0, left=0, right=1, top=1, bottom=0)
for i in np.arange(ny):
    for j in np.arange(nx):
        axs[i][j].set_xlim([-80,60])
        axs[i][j].step(v, cube[:,center_y+i,center_x+j], color='cornflowerblue', linewidth=2.)
        axs[i][j].plot(v, reconstructed_cube[:,center_y+i,center_x+j], linestyle="-", linewidth=2., color="k")
        iid = np.where((id1 == center_y+i) & (id2 == center_x+j))[0]
        for k in range(len(iid)):            
            axs[i][j].plot(v, gaussian(np.arange(len(v)), amp[iid[k]], mean[iid[k]], sigma[iid[k]]), linewidth=2., color=plt.cm.inferno(pvalues[k]))
        if j == 0: axs[i][j].set_ylabel(r'T [k]')
        axs[i][j].set_xlabel(r'v [km s$^{-1}$]')
plt.savefig("plot/" + 'mosaic_spectra_all_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

# #Plot centroid velocity field CNM                                                                                                                                                                                                                                                  
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_axes([0.1,0.1,0.74,0.8])
# ax.set_xlabel(r"x", fontsize=18.)
# ax.set_ylabel(r"y", fontsize=18.)
# ax.axes.xaxis.set_ticklabels([])
# ax.axes.yaxis.set_ticklabels([])
# ax.tick_params(labelsize=16)
# img = ax.imshow(vfield[9], aspect='auto', **imkw_coolwarm)
# divider = make_axes_locatable(ax)
# cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
# fig.add_axes(cax)
# cbar = fig.colorbar(img, cax=cax, orientation="horizontal", extend='both')
# cbar.ax.tick_params(labelsize=16.)
# cbar.set_label(r"Velocity Centroid [km.s$^{-1}$]", fontsize=18.)
# plt.savefig("plot/" + 'CV_LOCAL_CNM_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

# #Plot centroid velocity field INTER                                                                                                                                                                                                                                                  
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_axes([0.1,0.1,0.74,0.8])
# ax.set_xlabel(r"x", fontsize=18.)
# ax.set_ylabel(r"y", fontsize=18.)
# ax.axes.xaxis.set_ticklabels([])
# ax.axes.yaxis.set_ticklabels([])
# ax.tick_params(labelsize=16)
# img = ax.imshow(CV_INTER_LOCAL_ROHSA, aspect='auto', **imkw_coolwarm)
# divider = make_axes_locatable(ax)
# cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
# fig.add_axes(cax)
# cbar = fig.colorbar(img, cax=cax, orientation="horizontal", extend='both')
# cbar.ax.tick_params(labelsize=16.)
# cbar.set_label(r"Velocity Centroid [km.s$^{-1}$]", fontsize=18.)
# plt.savefig("plot/" + 'CV_LOCAL_INTER_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot integrated column density field
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.78,0.8], projection=w)
ax.set_xlabel(r"Gal. lon.", fontsize=18.)
ax.set_ylabel(r"Gal. lat.", fontsize=18.)
img = ax.imshow(vfield[9], aspect='auto', vmin=np.min(vfield[9]), vmax=np.max(vfield[9]), **imkw_coolwarm)
colorbar_ax = fig.add_axes([0.89, 0.1, 0.02, 0.8])
cbar = fig.colorbar(img, cax=colorbar_ax, extend='both')
cbar.ax.tick_params(labelsize=14.) 
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=18.)
plt.savefig("plot/" + 'CV_CNM_NEP.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)
