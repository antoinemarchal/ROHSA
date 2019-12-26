import numpy as np
from astropy.io import fits
from astropy import constants as const
from astropy import units as u
from scipy import optimize
from scipy import ndimage
from scipy import signal

from ROHSApy import ROHSA
import marchalib as ml

import function as fct
reload(fct)

#Open simulated data
path="/mnt/raid-cita/amarchal/SIMUSED/data/"
fitsname = "intensity_adim.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
cube = hdu[0].data

#NHI 
shape = (128,128)
NHI = ml.fBmnd(shape, ml.Pkgen(3.6,0.,np.inf), seed=27, unit_length=1)
NHI = np.interp(NHI, (NHI.min(), NHI.max()), (3., 4.))

#Open colour correction polynomial coefficient
color = fits.open("/home/amarchal/ROHSA/data/col_cor_iras_hfi_DX9v2_poly.fits")[0].data

#Parameters
n_mbb=1
l0=(849.27041926*u.micron).cgs.value
degree=5
lambda_tau = 1.
lambda_beta = 0.1
lambda_td = 0.1
lambda_var_sig = 2.
lambda_butter = 0.

kernel = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]]) / 4.     
cexp = (const.h.cgs*const.c.cgs/const.k_B.cgs).value
freq = np.array([np.full((cube.shape[1],cube.shape[2]),i) for i in np.array([3000,857,545,353])])*u.GHz
wavelength = (const.c / freq).to(u.micron).cgs.value

kmat = ml.kgrid(cube[0])
butter = np.sqrt(fct.butterworth(kmat,1.,0.06,-2.2))

theta = np.zeros((3*n_mbb,cube.shape[1],cube.shape[2]))
lbounds = np.zeros((3*n_mbb,cube.shape[1],cube.shape[2]))
ubounds = np.zeros((3*n_mbb,cube.shape[1],cube.shape[2]))

#init bounds and theta
lbounds[0] = 0. ; ubounds[0] = 100.; theta[0] = 0.5
lbounds[1] = 1. ; ubounds[1] = 2.5 ; theta[1] = 1.
lbounds[2] = 8. ; ubounds[2] = 50. ; theta[2] = 14.
# lbounds[3] = 0.   ; ubounds[3] = 100.
# lbounds[4] = 1.   ; ubounds[4] = 2.5
# lbounds[5] = 8.   ; ubounds[5] = 50.

bounds = [(lbounds.ravel()[i], ubounds.ravel()[i]) for i in np.arange(len(lbounds.ravel()))]

#bounds m1 variable
bounds.append([0.,10.]) 

#Concatanate all variables
n_theta = 3*n_mbb*cube.shape[1]*cube.shape[2]
params = np.zeros(n_theta+1)

params[:n_theta] = theta.ravel()
params[-1] = 0.41

def f_g(pars, n_pars, n_mbb, wave, data, l0, color):
    params = np.reshape(pars[:n_pars], (3*n_mbb, data.shape[1], data.shape[2]))
    m1 = pars[-1]

    grad_lin = np.zeros(pars.shape)
    dF_over_dB = np.zeros((params.shape[0], data.shape[0], data.shape[1], data.shape[2]))
    conv = np.zeros((params.shape[0], data.shape[1], data.shape[2]))
    dR_over_dB = np.zeros((params.shape[0], data.shape[1], data.shape[2]))

    #model Dust only
    # model = fct.MBB_l(wave,params[0], params[1], params[2],l0)
    # # dust = fct.MBB_l_adim(wave,params[3], params[4], params[5]*u.K,l0)

    # dF_over_dB[0] = fct.d_MBB_l_d_tau(wave,params[0],params[1],params[2],l0)
    # dF_over_dB[1] = fct.d_MBB_l_d_beta(wave,params[0],params[1],params[2],l0)
    # dF_over_dB[2] = fct.d_MBB_l_d_td(wave,params[0],params[1],params[2],l0)
    
    #Fast
    #precompute quantities
    tau = params[0]
    beta = params[1]
    td = params[2]
    l = wave

    pl = 1./(np.exp(cexp/l/td) - 1.)
    l0_l = l0/l
    a = l0_l**beta
    b = cexp / l
    b_td = b / td
    exp_b_td = np.exp(b_td)
    apl = a * pl

    #Attache aux donnees
    model =  tau * apl
    F = model - data

    dF_over_dB[0] = apl
    dF_over_dB[1] = tau * np.log(l0_l) * apl
    dF_over_dB[2] = tau * a * b * exp_b_td / (td**2. * (exp_b_td - 1.)**2.)

    d_F_times_F = np.sum(dF_over_dB * F,axis=1)

    #Smoothness Laplacian filtering
    conv[0] = lambda_tau * ndimage.convolve(tau, kernel, mode='reflect')
    conv[1] = lambda_beta * ndimage.convolve(beta, kernel, mode='reflect')
    conv[2] = lambda_td * ndimage.convolve(td, kernel, mode='reflect')
    
    dR_over_dB[0] = lambda_tau * ndimage.convolve(conv[0], kernel, mode='reflect')
    dR_over_dB[1] = lambda_beta * ndimage.convolve(conv[1], kernel, mode='reflect')
    dR_over_dB[2] = lambda_td * ndimage.convolve(conv[2], kernel, mode='reflect')

    #Variance / NHI correlation
    R_var_sig = lambda_var_sig * np.sum(((tau/NHI) - m1)**2.)
    grad_R_var_sig = lambda_var_sig * (tau/NHI - m1) / NHI

    #Butterworth filtering tau
    tf_tau = np.fft.fftshift(np.fft.fft2(tau))
    tf_tau_filtered = lambda_butter * np.fft.ifftn(np.fft.ifftshift((np.sqrt(butter) * np.abs(tf_tau))**2.)).real
    d_tf_tau_filtered_d_tau = lambda_butter * np.fft.ifftn(np.fft.ifftshift((np.conj(np.sqrt(butter)) * np.sqrt(butter) * tf_tau))).real

    #L2 norm / energy
    J = np.sum((model - data)**2)
    R = np.sum(conv**2)
    Q = np.sum(tf_tau_filtered)

    L = J + R + Q + R_var_sig
    grad = d_F_times_F + dR_over_dB 
    grad[0] += d_tf_tau_filtered_d_tau
    grad[0] += grad_R_var_sig

    grad_lin[:n_pars] = grad.ravel()
    grad_lin[-1] = - lambda_var_sig * np.sum(((tau/NHI) - m1))

    return 0.5*L, grad_lin
    

result = optimize.fmin_l_bfgs_b(f_g, params, 
                                args=(n_theta, n_mbb, wavelength, cube, 
                                      l0, color), bounds=bounds, 
                                approx_grad=False, disp=1, maxiter=800)
    
theta = np.reshape(result[0][:n_theta],(3*n_mbb, cube.shape[1], cube.shape[2]))


