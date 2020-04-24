import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from ROHSApy import ROHSA

path = "/data/amarchal/ROHSA_paper/"
fitsname = path + "data/synthetic_obs/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2.fits"
hdu = fits.open(fitsname)
hdr = hdu[0].header
cube = hdu[0].data

filename_parameters = path + "PARAMETERS/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_parameters_run_42.txt"
filename = path + "DAT/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2.dat"
fileout = path + "ROHSA/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_42.dat"
timeout = path + "ROHSA/Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_42_timestep.dat"
filename_noise = ""
n_gauss = 8
lambda_amp = 10000.                                                                                                                      
lambda_mu = 10000.
lambda_sig = 10000.                                                                               
lambda_var_amp = 0.                                                                                                          
lambda_var_mu = 0.                                                                                                        
lambda_var_sig = 1000. 
amp_fact_init = 0.66                                                                                                      
lb_sig_init = 0.001
ub_sig_init = 100.
lb_sig = 0.001
ub_sig = 100.
sig_init = 2.                                                                                                       
init_option = 'mean'                                                                                            
maxiter_init = 15000                                                                                              
maxiter = 800                                                                                   
m = 10                                                                                                
noise = ".false."                                                                                                              
regul = ".true."                                                                      
descent = ".true."                                                                                                                 
lstd = 1                                                                                                                        
ustd = 10                                                                                                                         
iprint = -1                                                                                                                     
iprint_init = -1                                                 
save_grid = ".true."  

core = ROHSA(cube)            
# core.cube2dat(filename=filename)
core.gen_parameters(filename_parameters=filename_parameters,
                    filename=filename, 
                    fileout=fileout,  
                    timeout=timeout,  
                    filename_noise = filename_noise,
                    n_gauss=n_gauss,
                    lambda_amp=lambda_amp,
                    lambda_mu=lambda_mu,
                    lambda_sig=lambda_sig,
                    lambda_var_sig=lambda_var_sig,
                    maxiter=maxiter,
                    amp_fact_init=amp_fact_init,
                    lb_sig=lb_sig,
                    ub_sig=ub_sig,
                    lb_sig_init=lb_sig_init,
                    ub_sig_init=ub_sig_init,
                    sig_init=sig_init,
                    noise=noise,
                    descent=descent,
                    lstd=lstd,
                    ustd=ustd,
                    iprint=iprint,
                    iprint_init=iprint_init,
                    save_grid=save_grid)
