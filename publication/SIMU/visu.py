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
cm_inferno.set_under(color='black')
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
data = np.genfromtxt("/data/amarchal/ROHSA_paper/ROHSA/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_34.dat") # 23 now 34 with new ROHSA

cube = hdu[0].data
cube[np.where(cube != cube)] = 0.

dim_x = cube.shape[2]
dim_y = cube.shape[1]
dim_v = cube.shape[0]

CDELT = hdr['CDELT3'] #*1.e-3  #km.s-1
CRVAL = hdr['CRVAL3'] #*1.e-3  #km.s-1
CRPIX = hdr['CRPIX3']

reso = np.abs(CDELT)

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

cube_INTER_ROHSA = np.zeros(cube.shape)
for i in np.arange(cube.shape[1]):
    for j in np.arange(cube.shape[2]):
        cube_INTER_ROHSA[:,i,j] = gaussian(np.arange(dim_v), params[0,i,j], vel2mean(CRVAL, CDELT, CRPIX, [params[1,i,j]])[0], params[2,i,j]/reso) + gaussian(np.arange(dim_v), params[3,i,j], vel2mean(CRVAL, CDELT, CRPIX, [params[4,i,j]])[0], params[5,i,j]/reso)

CV_INTER_ROHSA = np.tensordot(v, cube_INTER_ROHSA, axes=([0],[0])) / np.sum(cube_INTER_ROHSA, axis=0)

mean_v = [np.around(np.mean(vfields), decimals=1) for vfields in vfield]
mean_sig = [np.around(np.mean(sigfields), decimals=1) for sigfields in sigfield]

print mean_v
print mean_sig

print "max WNM", np.max(sigfield[0])
print "max LNM", np.max([np.max(sigfield[2]), np.max(sigfield[5])])
print "max CNM", np.max([np.max(sigfield[1]), np.max(sigfield[3]), np.max(sigfield[4])])

stop

# PLOT HIST SIGMA  / AMP                                                                                                                                                                                                                    
fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.tick_params(labelsize=14)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
bins = np.linspace(np.min(sigma), np.max(sigma), 400)
ax.hist(sigma, weights=np.sqrt(2.*np.pi)*amp*sigma*C/(np.sum(field)*1.e18), bins=bins, log=False, histtype='step', color='k', normed=False, 
        label=r'$\sigma weighted \, by \, A$', linewidth=2.5)
ax.set_xlim([0., 10])
# ax.set_ylim([0., 4.])  
ax.set_xlabel(r'$\sigma$ [km s$^{-1}$]', fontsize = 14)                                                                                  
ax.set_ylabel(r'PDF [fraction of total emission]', fontsize = 14)          
plt.savefig('plot/PDF_sigma_over_A.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02) 

#Good one
NHI_WNM = np.sqrt(2.*np.pi) * ampfield[0] * sigfield[0] * C / 1.e18
NHI_LNM = np.sqrt(2.*np.pi) * ampfield[1]*sigfield[1] * C / 1.e18
NHI_CNM = np.sqrt(2.*np.pi) * (ampfield[2]*sigfield[2] + ampfield[3]*sigfield[3] + ampfield[4]*sigfield[4] 
                               + ampfield[5]*sigfield[5] + ampfield[6]*sigfield[6] + ampfield[7]*sigfield[7]) * C / 1.e18
NHI_INTER = NHI_LNM + NHI_WNM

np.save("NHI_WNM_ROHSA.npy", NHI_WNM)
np.save("NHI_LNM_ROHSA.npy", NHI_LNM)
np.save("NHI_CNM_ROHSA.npy", NHI_CNM)

NHI_TOT = NHI_CNM + NHI_LNM + NHI_WNM

cube_CNM_simu = fits.open("/data/amarchal/Saury2014/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_0_Tmax_500_ROHSA_noise_0.0_K_beam_0_2_2.fits")[0].data
cube_LNM_simu = fits.open("/data/amarchal/Saury2014/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_500_Tmax_5000_ROHSA_noise_0.0_K_beam_0_2_2.fits")[0].data
cube_WNM_simu = fits.open("/data/amarchal/Saury2014/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_5000_Tmax_inf_ROHSA_noise_0.0_K_beam_0_2_2.fits")[0].data
cube_INTER_simu = fits.open("/data/amarchal/Saury2014/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_500_Tmax_inf_ROHSA_noise_0.0_K_beam_0_2_2.fits")[0].data

NHI_CNM_simu = np.sum(cube_CNM_simu,axis=0)
NHI_LNM_simu = np.sum(cube_LNM_simu,axis=0)
NHI_WNM_simu = np.sum(cube_WNM_simu,axis=0)
NHI_INTER_simu = np.sum(cube_INTER_simu,axis=0)

NHI_TOT_simu = NHI_CNM_simu + NHI_LNM_simu + NHI_WNM_simu

NHI_TOT_simu = NHI_TOT_simu * reso * C / 1.e18
NHI_CNM_simu = NHI_CNM_simu * reso * C / 1.e18
NHI_LNM_simu = NHI_LNM_simu * reso * C / 1.e18
NHI_WNM_simu = NHI_WNM_simu * reso * C / 1.e18
NHI_INTER_simu = NHI_INTER_simu * reso * C / 1.e18

# NHI_CNM_simu = np.load("NHI_CNM_paper_AM.npy")
# NHI_LNM_simu = np.load("NHI_LNM_paper_AM.npy")
# NHI_WNM_simu = np.load("NHI_WNM_paper_AM.npy")
# NHI_TOT_simu = NHI_CNM_simu + NHI_LNM_simu + NHI_WNM_simu

er_NHI =  np.abs((np.sum(NHI_TOT_simu) - np.sum(NHI_TOT)) / np.sum(NHI_TOT_simu)) * 100.

original_map = np.sum(cube, axis=0) * reso * C / 1.e18
reconstructed_map = np.sum(reconstructed_cube, axis=0) * reso * C / 1.e18

mu3 = np.abs(moment((reconstructed_cube - cube).ravel(), moment=3, nan_policy='propagate'))
print "mu3 = ", mu3

#SPS1D
stat_cnm = tb.PowerS(NHI_CNM)
stat_inter = tb.PowerS(NHI_INTER)

stat_cnm_simu = tb.PowerS(NHI_CNM_simu)
stat_inter_simu = tb.PowerS(NHI_INTER_simu)

ks = stat_inter_simu.get_ks(unit_length=1)

sps1d_cnm = stat_cnm.sps1d(return_log=False)
sps1d_inter = stat_inter.sps1d(return_log=False)

sps1d_cnm_simu = stat_cnm_simu.sps1d(return_log=False)
sps1d_inter_simu = stat_inter_simu.sps1d(return_log=False)

ksup = np.where(np.array(ks) > 0.4)[0][0]

# #SPS1D CENTROID VELOCITY
# CV_WNM_simu = fits.open("/data/amarchal/Saury2014/synthetic_obs/CV_WNM_paper.fits")[0].data
CV_INTER_simu = np.tensordot(v, cube_INTER_simu, axes=([0],[0])) / np.sum(cube_INTER_simu, axis=0)

stat_CV_inter = tb.PowerS(CV_INTER_ROHSA)
stat_CV_inter_simu = tb.PowerS(CV_INTER_simu)

sps1d_CV_inter = stat_CV_inter.sps1d(return_log=False)
sps1d_CV_inter_simu  = stat_CV_inter_simu.sps1d(return_log=False)

stop

def plot_spect(x,y,velocity,cube,id1,id2,amp,mean,sigma,reso):
    clr = ['y', 'm', 'g', 'b', 'orange', 'cyan', 'y', 'm', 'g', 'b', 'orange', 'cyan',
           'y', 'm', 'g', 'b', 'orange', 'cyan', 'y', 'm', 'g', 'b', 'orange', 'cyan']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.step(velocity, cube[:,x,y], color='cornflowerblue')
    ax.plot(velocity, reconstructed_cube[:,x,y], color='r')
    iid = np.where((id1 == x) & (id2 == y))[0]
    for i in range(len(iid)):
        ax.plot(velocity, gaussian(np.arange(len(v)), amp[iid[i]], mean[iid[i]], sigma[iid[i]]/reso), color=clr[i])
    for i in range(len(iid)):
        print amp[iid[i]], vmean[iid[i]], sigma[iid[i]]

#Plot integrated column density field
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.74,0.8])
ax.set_xlabel(r"x", fontsize=18.)
ax.set_ylabel(r"y", fontsize=18.)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
# ax.axis('off')
ax.tick_params(labelsize=16)
img = ax.imshow(original_map, aspect='auto', vmin=np.min(original_map), vmax=np.max(original_map), **imkw_inferno)
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
fig.add_axes(cax)
cbar = fig.colorbar(img, cax=cax, orientation="horizontal", extend='both')
cbar.ax.tick_params(labelsize=16.) 
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)
plt.savefig("plot/" + 'NHI_TOT.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot mosaic field
lh = 4; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,23.))
fig.subplots_adjust(top=1., bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
for i in np.arange(lh):
    for j in np.arange(lw):
        im = axs[i][j].imshow(field[k], **imkw_inferno)
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        # axs[i][j].set_title(r"$\overline{\mu}$ = " + str(np.around(np.mean(vfield[k]), decimals=1)) + 
        #                     ", $\overline{\sigma}$ = " + str(np.around(np.mean(sigfield[k]), decimals=1)), fontsize=10.)
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='4%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%d', extend='both')
        cbar.ax.tick_params(labelsize=16.) 
        if i == lh-1 : cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)
        k += 1
plt.savefig('plot/mosaic_field.pdf', format='pdf')

#Plot mosaic vfield
# lh = 3; lw = 2
# fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,17.))
lh = 4; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,23.))
fig.subplots_adjust(top=1., bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
mincbar = [None, None, None, None, None, None, None, None]
maxcbar = [None, None, None, None, None, None, None, None]
for i in np.arange(lh):
    for j in np.arange(lw):
        # im = axs[i][j].imshow(field[k], origin="lower", interpolation=None, cmap="gist_gray")
        im1 = axs[i][j].imshow(vfield[k], vmin=mincbar[k] , vmax=maxcbar[k], **imkw_coolwarm)
        # im2 = axs[i][j].contour(field[k], colors='k', linestyles='-', linewidths=0.2)
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='4%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
        cbar.ax.tick_params(labelsize=16.) 
        if i == lh-1 : cbar.set_label(r"v [km s$^{-1}$]", fontsize=18.)
        k += 1
plt.savefig('plot/mosaic_vfield.pdf', format='pdf')

#Plot mosaic sigfield
# lh = 3; lw = 2
# fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,17.))
lh = 4; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,23.))
fig.subplots_adjust(top=1., bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
k = 0
# mincbar = [None, 1.1, 5., 1.6, None, 4.65]
# maxcbar = [None, 2.3, 5.6, None, None, None]
mincbar = [None, None, None, None, None, None, None, None]
maxcbar = [None, None, None, None, None, None, None, None]
for i in np.arange(lh):
    for j in np.arange(lw):
        # im = axs[i][j].imshow(field[k], origin="lower", interpolation=None, cmap="gist_gray")
        im1 = axs[i][j].imshow(sigfield[k], vmin=mincbar[k] , vmax=maxcbar[k], **imkw_cubehelix)
        # im2 = axs[i][j].contour(field[k], colors='k', linestyles='-', linewidths=0.2)
        if j == 0: axs[i][j].set_ylabel(r'y')
        axs[i][j].set_xlabel(r'x')
        axs[i][j].axes.xaxis.set_ticklabels([])
        axs[i][j].axis('off')
        divider = make_axes_locatable(axs[i][j])
        cax = divider.append_axes('bottom', size='4%', pad=0.05)
        cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
        cbar.ax.tick_params(labelsize=16.) 
        if i == lh-1 : cbar.set_label(r"$\sigma_v}$ [km s$^{-1}$]", fontsize=18.)
        k += 1
plt.savefig('plot/mosaic_sigfield.pdf', format='pdf')

#Plot mosaic spectra                                                                                                                    
pvalues = np.logspace(-1, 0, 10)
pmin = pvalues[0]
pmax = pvalues[-1]

def norm(pval):
    return (pval - pmin) / float(pmax - pmin)

ny = 4; nx = 4                                                                                                                                          
center_y = 36; center_x = 42
cb = "magenta"
cw = "crimson"
fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(cm2inch((18.,14.))))
fig.subplots_adjust(hspace=0, wspace=0, left=0, right=1, top=1, bottom=0)
for i in np.arange(ny):
    for j in np.arange(nx):
        # axs[i][j].set_yticks(np.arange(0., 18., 5.))
        # axs[i][j].set_xticks(np.arange(-10, 20., 10.))
        # axs[i][j].set_ylim([-1,20])
        axs[i][j].set_xlim([-25,25])
        axs[i][j].tick_params(labelsize=8)
        axs[i][j].step(v, cube[:,center_y+i,center_x+j], color='cornflowerblue', linewidth=2.)
        axs[i][j].plot(v, reconstructed_cube[:,center_y+i,center_x+j], linestyle="-", linewidth=2., color="k")
        iid = np.where((id1 == center_y+i) & (id2 == center_x+j))[0]
        for k in range(len(iid)):            
            axs[i][j].plot(v, gaussian(np.arange(len(v)), amp[iid[k]], mean[iid[k]], np.sqrt((sigma[iid[k]]/reso)**2. + 1.**2)), linewidth=2., color=plt.cm.inferno(norm(pvalues[k])))
        if j == 0: axs[i][j].set_ylabel(r'T [k]', fontsize=8)
        axs[i][j].set_xlabel(r'v [km s$^{-1}$]', fontsize=8)
plt.savefig("plot/" + 'mosaic_spectra_all.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#HEATMAP
x_bins = np.linspace(np.min(vmean), np.max(vmean), 800)
y_bins = np.logspace(np.log(0.1), np.log(np.max(sigma)), 800)
H, xedges, yedges = np.histogram2d(vmean, sigma, weights=np.sqrt(2.*np.pi)*amp*sigma, bins=[x_bins, y_bins])
H = np.ma.masked_invalid(np.atleast_2d(H))
fig = plt.figure(figsize=(16.,8.))
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_ylim([1., 20.])
ax.set_xlim([np.min(vmean),np.max(vmean)])
ax.set_xlabel(r'$v_r [km.s^{-1}]$',  fontsize = 16)
ax.set_ylabel(r'$\sigma [km.s^{-1}]$',  fontsize = 16)
ax.yaxis.set_major_formatter(ScalarFormatter()) 
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='2%', pad=0.05)
im = ax.pcolormesh(xedges, yedges, np.log(H.T), vmin=0, vmax=np.max(np.log(H.T)), cmap=cm_inferno)
cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.set_label(r'$log(W_{HI}) [K.km.s^{-1}]$', fontsize = 16)
plt.savefig("plot/" + 'heatmap.png', format='png', bbox_inches='tight', pad_inches=0.02)

#Plot mosaic field
lh = 3; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,17.))
fig.subplots_adjust(top=1., bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
# fig.subplots_adjust(top=1.015, bottom=0.01, left=0., right=1., hspace=0., wspace=0.02)

im = axs[0][0].imshow(NHI_WNM, vmin=np.min(NHI_WNM_simu), vmax=np.max(NHI_WNM_simu), **imkw_inferno)
axs[0][0].axes.xaxis.set_ticklabels([])
axs[0][0].axis('off')
divider = make_axes_locatable(axs[0][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 

im = axs[1][0].imshow(NHI_LNM, vmin=np.min(NHI_LNM_simu), vmax=np.max(NHI_LNM_simu), **imkw_inferno)
axs[1][0].axes.xaxis.set_ticklabels([])
axs[1][0].axis('off')
divider = make_axes_locatable(axs[1][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 

im = axs[2][0].imshow(NHI_CNM, vmin=np.min(NHI_CNM_simu), vmax=np.max(NHI_CNM_simu), **imkw_inferno)
axs[2][0].axes.xaxis.set_ticklabels([])
axs[2][0].axis('off')
divider = make_axes_locatable(axs[2][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

im = axs[0][1].imshow(NHI_WNM_simu, **imkw_inferno)
axs[0][1].axes.xaxis.set_ticklabels([])
axs[0][1].axis('off')
divider = make_axes_locatable(axs[0][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 

im = axs[1][1].imshow(NHI_LNM_simu, **imkw_inferno)
axs[1][1].axes.xaxis.set_ticklabels([])
axs[1][1].axis('off')
divider = make_axes_locatable(axs[1][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 

im = axs[2][1].imshow(NHI_CNM_simu, **imkw_inferno)
axs[2][1].axes.xaxis.set_ticklabels([])
axs[2][1].axis('off')
divider = make_axes_locatable(axs[2][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.tick_params(labelsize=16) 
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

plt.savefig('plot/mosaic_field_comparison.pdf', format='pdf')

stop

# fig = plt.figure(figsize=(cm2inch((18.,18.))))
# ax = fig.add_subplot(111)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.xaxis.set_major_formatter(ScalarFormatter())
# ax.yaxis.set_major_formatter(ScalarFormatter())
# bins_wnm = np.logspace(np.log10(0.01), np.log10(np.nanmax(NHI_WNM.ravel())), 200)
# ax.hist(NHI_WNM.ravel(), bins=bins_wnm, log=True, histtype='step', color='r', normed=False, label='', linewidth=1.5)
# ax.hist(NHI_WNM_simu.ravel(), bins=bins_wnm, log=True, histtype='step', color='o', normed=False, label='', linewidth=1.5)

# bins_lnm = np.logspace(np.log10(0.01), np.log10(np.nanmax(NHI_LNM.ravel())), 200)
# ax.hist(NHI_LNM.ravel(), bins=bins_lnm, log=True, histtype='step', color='g', normed=False, label='', linewidth=1.5)
# ax.hist(NHI_LNM_simu.ravel(), bins=bins_lnm, log=True, histtype='step', color='k', normed=False, label='', linewidth=1.5)

# bins_cnm = np.logspace(np.log10(0.01), np.log10(np.nanmax(NHI_CNM.ravel())), 200)
# ax.hist(NHI_CNM.ravel(), bins=bins_cnm, log=True, histtype='step', color='g', normed=False, label='', linewidth=1.5)
# ax.hist(NHI_CNM_simu.ravel(), bins=bins_cnm, log=True, histtype='step', color='k', normed=False, label='', linewidth=1.5)
# # ax.set_xlim([0.01, np.max(NHI_theo_list)])
# # ax.set_ylim([1., 1.e5])
# for axis in [ax.xaxis]:
#     formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
#     axis.set_major_formatter(formatter)
# plt.setp(ax.get_yticklabels()[1], visible=False)
# plt.setp(ax.get_xticklabels()[1], visible=False)
# ax.set_xlabel(r'$NHI$ [1.83 10$^{18}$ cm$^{-2}$]', fontsize=14)
# ax.set_ylabel(r'Normalized $PDF$', fontsize=14)
# plt.savefig('plot/PDF_NHI.pdf', format='pdf')

#SPS1D NHI
fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
# ax.plot(ks[:ksup], sps1d_wnm[:ksup], linestyle='-', color="red", linewidth=3.5)
# ax.plot(ks[:ksup], sps1d_wnm_simu[:ksup], '.', color='orange', markersize=8)
# ax.plot(ks[:ksup], sps1d_lnm[:ksup], linestyle='-', color="k", linewidth=3.5)
# ax.plot(ks[:ksup], sps1d_lnm_simu[:ksup], '.k', markersize=8) 
ax.plot(ks[:ksup], sps1d_inter[:ksup], linestyle='-', color="red", linewidth=3.5)
ax.plot(ks[:ksup], sps1d_inter_simu[:ksup], '.', color="orange", markersize=8) 
ax.plot(ks[:ksup], sps1d_cnm[:ksup], linestyle='-', color="blue", linewidth=3.5)
ax.plot(ks[:ksup], sps1d_cnm_simu[:ksup], '.', color="cyan", markersize=8) 
# ax.plot(ks[:ksup], sps1d_tot[:ksup], linestyle='-', color="black", linewidth=3.5)
# ax.plot(ks[:ksup], sps1d_tot_simu[:ksup], '.', color="black", markersize=8) 
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim([0.005, 0.5])
ax.set_ylim([1.e4, 1.e12])
for axis in [ax.xaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)
plt.setp(ax.get_yticklabels()[1], visible=False)
plt.setp(ax.get_xticklabels()[1], visible=False)

ax.set_xlabel(r'k [pixel$^{-1}$]',  fontsize = 16)
ax.set_ylabel(r'P(k) [Arbitrary unit]',  fontsize = 16)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.savefig("plot/sps1d_NHI.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)

#ER SPS1D NHI
fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.set_ylim([-1., 1.])
ax.plot(ks[:ksup], (np.array(sps1d_wnm_simu[:ksup])-np.array(sps1d_wnm[:ksup]))/np.array(sps1d_wnm_simu[:ksup]), linestyle='-', color="red", linewidth=3.5)
ax.plot(ks[:ksup], (np.array(sps1d_lnm_simu[:ksup])-np.array(sps1d_lnm[:ksup]))/np.array(sps1d_lnm_simu[:ksup]), linestyle='-', color="black", linewidth=3.5)
ax.plot(ks[:ksup], (np.array(sps1d_cnm_simu[:ksup])-np.array(sps1d_cnm[:ksup]))/np.array(sps1d_cnm_simu[:ksup]), linestyle='-', color="blue", linewidth=3.5) 
ax.set_xlabel(r'k [pixel$^{-1}$]',  fontsize = 16)
ax.set_ylabel(r'P$_{SIMU}$(k) - P$_{ROHSA}$(k) / P$_{SIMU}$(k)',  fontsize = 16)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
plt.savefig("plot/sps1d_NHI_rapport.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)

#Plot centroid velocity field INTER simu
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.74,0.8])
ax.set_xlabel(r"x", fontsize=18.)
ax.set_ylabel(r"y", fontsize=18.)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(labelsize=16)
img = ax.imshow(CV_INTER_simu, aspect='auto', **imkw_coolwarm)
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
fig.add_axes(cax)
cbar = fig.colorbar(img, cax=cax, orientation="horizontal", extend='both')
cbar.ax.tick_params(labelsize=16.) 
cbar.set_label(r"Velocity Centroid [km s$^{-1}$] SIMU", fontsize=18.)
plt.savefig("plot/" + 'CV_INTER_simu.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

#Plot centroid velocity field INTER ROHSA
fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.1,0.1,0.74,0.8])
ax.set_xlabel(r"x", fontsize=18.)
ax.set_ylabel(r"y", fontsize=18.)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(labelsize=16)
img = ax.imshow(CV_INTER_ROHSA, aspect='auto', vmin=np.min(CV_INTER_simu), vmax=np.max(CV_INTER_simu), **imkw_coolwarm)
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
fig.add_axes(cax)
cbar = fig.colorbar(img, cax=cax, orientation="horizontal", extend='both')
cbar.ax.tick_params(labelsize=16.) 
cbar.set_label(r"Velocity Centroid [km s$^{-1}$] ROHSA", fontsize=18.)
plt.savefig("plot/" + 'CV_INTER_ROHSA.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)

lh = 1; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10.,5.8))
fig.subplots_adjust(top=1.1, bottom=0.0, left=0.01, right=0.99, hspace=0.01, wspace=0.02)
mincbar = [np.min(CV_INTER_simu), None]
maxcbar = [np.max(CV_INTER_simu), None]
im = axs[0].imshow(CV_INTER_ROHSA, vmin=mincbar[0] , vmax=maxcbar[0], **imkw_coolwarm)
axs[0].set_ylabel(r'y')
axs[0].set_xlabel(r'x')
axs[0].axes.xaxis.set_ticklabels([])
axs[0].axis('off')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.) 
cbar.set_label(r"Velocity Centroid [km s$^{-1}$] ROHSA", fontsize=16.)

im = axs[1].imshow(CV_INTER_simu, vmin=mincbar[1] , vmax=maxcbar[1], **imkw_coolwarm)
axs[1].set_ylabel(r'y')
axs[1].set_xlabel(r'x')
axs[1].axes.xaxis.set_ticklabels([])
axs[1].axis('off')
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.) 
cbar.set_label(r"Velocity Centroid [km s$^{-1}$] SIMU", fontsize=16.)
plt.savefig('plot/CV_INTER_long.pdf', format='pdf')

lh = 2; lw = 1
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(5,11.))
fig.subplots_adjust(top=1.05, bottom=0.04, left=0.01, right=0.99, hspace=0.0, wspace=0.02)
mincbar = [np.min(CV_INTER_simu), None]
maxcbar = [np.max(CV_INTER_simu), None]
im = axs[0].imshow(CV_INTER_ROHSA, vmin=mincbar[0] , vmax=maxcbar[0], **imkw_coolwarm)
axs[0].set_ylabel(r'y')
axs[0].set_xlabel(r'x')
axs[0].axes.xaxis.set_ticklabels([])
axs[0].axis('off')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.)

im = axs[1].imshow(CV_INTER_simu, vmin=mincbar[1] , vmax=maxcbar[1], **imkw_coolwarm)
axs[1].set_ylabel(r'y')
axs[1].set_xlabel(r'x')
axs[1].axes.xaxis.set_ticklabels([])
axs[1].axis('off')
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=16.)
plt.savefig('plot/CV_INTER_lat.pdf', format='pdf')

#__________________________________________________________________________________________________________
lh = 2; lw = 2
fig, axs = plt.subplots(lh, lw, sharex=True, sharey=True, figsize=(10,12.))
fig.subplots_adjust(top=1.02, bottom=0.03, left=0.01, right=0.99, hspace=0.01, wspace=0.02)

im = axs[0][0].imshow(NHI_INTER, vmin=np.min(NHI_INTER_simu) , vmax=np.max(NHI_INTER_simu), **imkw_inferno)
axs[0][0].set_ylabel(r'y')
axs[0][0].set_xlabel(r'x')
axs[0][0].axes.xaxis.set_ticklabels([])
axs[0][0].axis('off')
divider = make_axes_locatable(axs[0][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%d', extend='both')
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

im = axs[0][1].imshow(NHI_INTER_simu, **imkw_inferno)
axs[0][1].set_ylabel(r'y')
axs[0][1].set_xlabel(r'x')
axs[0][1].axes.xaxis.set_ticklabels([])
axs[0][1].axis('off')
divider = make_axes_locatable(axs[0][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%d', extend='both')
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"N$_{HI}$ / [10$^{18}$ cm$^{-2}$]", fontsize=18.)

im1 = axs[1][0].imshow(CV_INTER_ROHSA, vmin=np.min(CV_INTER_simu) , vmax=np.max(CV_INTER_simu), **imkw_coolwarm)
# im2 = axs[1][0].contour(NHI_INTER, colors='k', linestyles='-', linewidths=0.2)
axs[1][0].set_ylabel(r'y')
axs[1][0].set_xlabel(r'x')
axs[1][0].axes.xaxis.set_ticklabels([])
axs[1][0].axis('off')
divider = make_axes_locatable(axs[1][0])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=16.)

im1 = axs[1][1].imshow(CV_INTER_simu, **imkw_coolwarm)
# im2 = axs[1][1].contour(NHI_INTER_simu, colors='k', linestyles='-', linewidths=0.2)
axs[1][1].set_ylabel(r'y')
axs[1][1].set_xlabel(r'x')
axs[1][1].axes.xaxis.set_ticklabels([])
axs[1][1].axis('off')
divider = make_axes_locatable(axs[1][1])
cax = divider.append_axes('bottom', size='4%', pad=0.05)
cbar = fig.colorbar(im1, cax=cax, orientation='horizontal', format='%.1f', extend='both')
cbar.ax.tick_params(labelsize=16.)
cbar.set_label(r"Velocity Centroid [km s$^{-1}$]", fontsize=16.)
plt.savefig('plot/NHI_CV_INTER.pdf', format='pdf')

#SPS1D CV
fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(ks[:ksup], np.array(sps1d_CV_inter[:ksup]), '.', color="orange", markersize=8)
ax.plot(ks[:ksup], np.array(sps1d_CV_inter_simu[:ksup]), linestyle='-', color='red', linewidth=3.5)
ax.set_xlabel(r'k [pixel$^{-1}$]',  fontsize = 16)
ax.set_ylabel(r'P(k) [Arbitrary unit]',  fontsize = 16)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
plt.savefig("plot/sps1d_CV_INTER.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)

fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.set_ylim([-1., 1.])
ax.plot(ks[:ksup], (np.array(sps1d_CV_wnm_simu[:ksup]) - np.array(sps1d_CV_wnm[:ksup]))/np.array(sps1d_CV_wnm_simu[:ksup]), linestyle='-', color='r', linewidth=3.5)
ax.set_xlabel(r'k [pixel$^{-1}$]',  fontsize = 16)
ax.set_ylabel(r'P$_{SIMU}$(k) - P$_{ROHSA}$(k) / P$_{SIMU}$(k)',  fontsize = 16)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
ax.xaxis.set_tick_params(labelsize=14)
ax.yaxis.set_tick_params(labelsize=14)
plt.savefig("plot/sps1d_CV_rapport.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)

#OBJ FUNCTION LAST ITERATION
obj_f = np.loadtxt("iterate_Tb_reso_0.807km.s-1_noise_0.01K_1024_beam_2_256_gauss_run_13.dat")

fig = plt.figure(figsize=(cm2inch((18.,18.))))
ax = fig.add_subplot(111)
ax.plot(np.arange(len(obj_f[:,9])), np.log10(obj_f[:,9]), linewidth=3.5)
ax.set_xlabel(r'Iteration',  fontsize = 16)
ax.set_ylabel(r'log$_{10}$ J($\bf \theta, \bf b$)',  fontsize = 16)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)
plt.savefig("plot/obf_f.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)
