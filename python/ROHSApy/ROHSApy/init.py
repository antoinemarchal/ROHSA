import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from ROHSApy import ROHSA

path = "/data/amarchal/"
fitsname = path + "data/observation/data.fits"
hdu = fits.open(fitsname)
hdr = hdu[0].header
cube = hdu[0].data[0]
# hdu_rms = fits.open(path + "data/")
# rms = hdu_rms[0].data

filename_parameters = path + "PARAMETERS/data_parameters_run_0.txt"
filename = path + "DAT/data.dat"
fileout = path + "ROHSA/data_gauss_run_0.dat"
# filename_noise = path + "DAT/data_rms.dat"
n_gauss = 12
lambda_amp = 1000.                                                                                                                      
lambda_mu = 1000.
lambda_sig = 1000.                                                                              
lambda_var_amp = 0.                                                                                                          
lambda_var_mu = 0.                                                                                                        
lambda_var_sig = 1000. 
amp_fact_init = 0.6666                                                                                                      
sig_init = 5.                                                                                                       
lb_sig_init = 0.001                                                                                                       
ub_sig_init = 12.                                                                                                       
lb_sig = 0.001                                                                                                       
ub_sig = 100.                                                                                                       
init_option = 'mean'                                                                                            
maxiter_init = 15000                                                                                              
maxiter = 800                                                                                   
m = 10                                                                                                
noise = ".false."                                                                                                              
regul = ".true."                                                                      
descent = ".true."                                                                                                                 
lstd = 1                                                                                                                        
ustd = 40                                                                                                                         
iprint = -1                                                                                                                     
iprint_init = -1                                                 
save_grid = ".true."  

core = ROHSA(cube)            

# core.cube2dat(filename=filename)
# core.rms_map(rms, filename=filename_noise)
core.gen_parameters(filename_parameters=filename_parameters,
                    filename=filename, 
                    fileout=fileout,  
                    # filename_noise = filename_noise,
                    n_gauss=n_gauss,
                    amp_fact_init=amp_fact_init,
                    lambda_amp=lambda_amp,
                    lambda_mu=lambda_mu,
                    lambda_sig=lambda_sig,
                    lambda_var_sig=lambda_var_sig,
                    lambda_var_mu=lambda_var_mu,
                    maxiter=maxiter,
                    sig_init=sig_init,
                    lb_sig=lb_sig,
                    ub_sig=ub_sig,
                    lb_sig_init=lb_sig_init,
                    ub_sig_init=ub_sig_init,
                    noise=noise,
                    descent=descent,
                    lstd=lstd,
                    ustd=ustd,
                    iprint=iprint,
                    iprint_init=iprint_init,
                    save_grid=save_grid)

