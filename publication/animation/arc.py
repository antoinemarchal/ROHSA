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
import pylab 
import random
from matplotlib.patches import Arc

import turbulence as tb

# plt.ion()

cm_coolwarm = plt.get_cmap('coolwarm')
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

def ellipse(x0,y0,a,b,alpha,phi):
    r = a*b/np.sqrt((b*np.cos(phi))**2 + (a*np.sin(phi))**2)
    return [x0+r*np.cos(phi+alpha), y0+r*np.sin(phi+alpha)]

C = 1.82243e18 

#DATA
hdu = fits.open("/data/amarchal/GHIGLS/data/GHIGLS_NEP_Tb.fits") 
hdr = hdu[0].header

cube = hdu[0].data[0]
cube[np.where(cube != cube)] = 0.

cube = cube[200:500,:,:]

dim_x = cube.shape[2]
dim_y = cube.shape[1]
dim_v = cube.shape[0]

CDELT = hdr['CDELT3'] *1.e-3  #km.s-1                                                                                                      
CRVAL = hdr['CRVAL3'] *1.e-3  #km.s-1                                                                                                             
CRPIX = hdr['CRPIX3'] - 200

reso = np.abs(CDELT)

vel = mean2vel(CRVAL, CDELT, CRPIX, np.arange(cube.shape[0]))

#FULL FIELD
data = np.genfromtxt("/data/amarchal/GHIGLS/data/DAT/ROHSA/GHIGLS_NEP_Tb_gauss_run_15.dat")

id1 = data[:, 0]
id2 = data[:, 1]
sigma = data[:, 4]
amp = data[:, 2]
mean = data[:, 3] - 1

vmean = mean2vel(CRVAL, CDELT, CRPIX, mean)

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

ys = 100
xs = np.arange(200)+10

# Figure setup
minvel_idx = np.where(np.array(vel[::-1]) > -80)[0][0]
maxvel_idx = np.where(np.array(vel[::-1]) < 20)[0][::-1][0]
subvel = vel[::-1][minvel_idx:maxvel_idx]

ang = np.interp(subvel, (np.min(subvel), np.max(subvel)), (20, 340))
color = ["b", "r", "k", "c", "g", "m", "y", "orange", "crimson", "k", "r", "b", "r", "k", "c", "g", "m", "y", "orange", "crimson", "k", "r"]

#Plot mosaic spectra random walk 2D 
for k in np.arange(len(xs)):
    ny = 4; nx = 4                                                                                                                                         
    center_y = ys; center_x = xs[k]
    cb = "magenta"
    cw = "crimson"

    fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(cm2inch((18.,14.))))
    fig.subplots_adjust(hspace=0, wspace=0, left=0, right=1, top=1, bottom=0)
    for i in np.arange(ny):
        for j in np.arange(nx):
            axs[i][j].set_xlim([-2,2])
            axs[i][j].set_ylim([-2,2])
            vmin = np.min([v[center_y+i,center_x+j] for v in vfield])
            vmax = np.max([v[center_y+i,center_x+j] for v in vfield])
    
            ampmax = np.max([amp[center_y+i,center_x+j] for amp in ampfield])

            for r in np.arange(n_gauss-1)+1:
                amplitude = ampfield[r][ys,xs[k]]
                velocity = vfield[r][ys,xs[k]]
                sig = 8.*sigfield[r][ys,xs[k]]
                
                idx = (np.abs(np.array(subvel)-velocity)).argmin()    
                ang_velocity = ang[idx]
                
                if amplitude == 0 : continue
                axs[i][j].add_patch(Arc((0., 0.), amplitude, amplitude, theta1=int(ang_velocity-sig), theta2=int(ang_velocity+sig), edgecolor=color[r-1], lw=1.5))
                x1,y1 = ellipse(0, 0, amplitude/2., amplitude/2., np.deg2rad(0), np.deg2rad(int(ang_velocity-sig)))
                x2,y2 = ellipse(0, 0, amplitude/2., amplitude/2., np.deg2rad(0), np.deg2rad(int(ang_velocity+sig)))
                axs[i][j].plot([0,x1],[0,y1], color=color[r-1],lw=1.5)
                axs[i][j].plot([0,x2],[0,y2], color=color[r-1],lw=1.5)
                axs[i][j].axes.get_xaxis().set_visible(False)
                axs[i][j].axes.get_yaxis().set_visible(False)
                axs[i][j].set_frame_on(False)

    plt.savefig("plot/mosaic_spectra_arc/" + 'mosaic_spectra_all_arc_NEP_{}.png'.format(k), format='png', bbox_inches='tight', pad_inches=0.02)

stop

for k in np.arange(len(xs)):
    fig_width, fig_height = 10, 10
    fig = plt.figure(figsize=(fig_width, fig_height), frameon=False)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], aspect='equal')
    ax.set_axis_off()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    
    vmin = np.min([v[ys,xs[k]] for v in vfield])
    vmax = np.max([v[ys,xs[k]] for v in vfield])
    
    ampmax = np.max([amp[ys,xs[k]] for amp in ampfield])

    for i in np.arange(n_gauss-1)+1:
        amplitude = ampfield[i][ys,xs[k]]
        velocity = vfield[i][ys,xs[k]]
        sig = 8.*sigfield[i][ys,xs[k]]
        
        idx = (np.abs(np.array(subvel)-velocity)).argmin()    
        ang_velocity = ang[idx]
        
        if amplitude == 0 : continue
        ax.add_patch(Arc((0., 0.), amplitude, amplitude, theta1=int(ang_velocity-sig), theta2=int(ang_velocity+sig), edgecolor=color[i-1], lw=1.5))
        x1,y1 = ellipse(0, 0, amplitude/2., amplitude/2., np.deg2rad(0), np.deg2rad(int(ang_velocity-sig)))
        x2,y2 = ellipse(0, 0, amplitude/2., amplitude/2., np.deg2rad(0), np.deg2rad(int(ang_velocity+sig)))
        ax.plot([0,x1],[0,y1], color=color[i-1],lw=1.5)
        ax.plot([0,x2],[0,y2], color=color[i-1],lw=1.5)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_frame_on(False)
    plt.savefig("plot/arc/" + 'arc_{}.png'.format(k), format='png', bbox_inches='tight', pad_inches=0.02)
