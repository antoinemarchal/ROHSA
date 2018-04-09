import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import matplotlib.colors as colors 
from astropy.io import fits
from astropy import units
from astropy import constants as const
import scipy.integrate as integrate

plt.ion()
cm = plt.get_cmap('inferno')
cm.set_bad(color='black')
imkw = dict(origin='lower', interpolation='none', cmap=cm)


path_simu = '/data/glx-calcul3/data1/amarchal/SIMULATIONS/Saury2014/n02_pw05_vs07/synthetic_obs/'

filename_cube = 'Tb_reso_0.2km.s-1_noise_0.6K_64.fits'
hdu_list_data = fits.open(path_simu + filename_cube)
cube = hdu_list_data[0].data

nx = cube.shape[2]
ny = cube.shape[1]
nv = cube.shape[0]
with open('cube.dat','w+') as f:
    f.write('{:d}\t{:d}\t{:d}\n'.format(nv, ny, nx))
    for i in range(nv):
        for j in range(ny):
            for k in range(nx):
                line = '{:d}\t{:d}\t{:d}\t{:5.8f}\n'.format(i, j, k, cube[i,j,k])
                f.write(line)
                
