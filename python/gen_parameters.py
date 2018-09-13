import os
import numpy as np

filename = ''
fileout = ''
filename_noise = ''
n_gauss = 6            # number of gaussian to fit
lambda_amp = 1.      # lambda for amplitude parameter
lambda_mu = 1.       # lamnda for mean position parameter
lambda_sig = 1.       # lambda for dispersion parameter
lambda_var_amp = 0.  # lambda for dispersion amplitude parameter
lambda_var_mu = 0.   # lambda for dispersion mean position parameter   
lambda_var_sig = 1.  # lambda for dispersion variance parameter
amp_fact_init = 0.66 # times max amplitude of additional Gaussian
sig_init = 5.        # dispersion of additional Gaussian
init_option = 'mean'   # init ROHSA with the mean or the std spectrum
maxiter_init = 15000   # max iteration for L-BFGS-B alogorithm (init mean spectrum)
maxiter = 800          # max iteration for L-BFGS-B alogorithm
m = 10                 # number of corrections used in the limited memory matrix by LBFGS-B
noise = ".false."        # if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
regul = ".true."         # if true --> activate regulation
descent = ".true."       # if true --> activate hierarchical descent to initiate the optimization
lstd = 1               # lower bound to compute the standard deviation map of the cube (if noise .eq. false)
ustd = 20              # upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
iprint = -1            # print option 
iprint_init = -1       # print option init

path = "./INPUT/"
if not os.path.exists(path):
    os.makedirs(path)

k = 0

input_file = open(path + "input_parameters_run_{}".format(k), 'w')
input_file.write("&user_parameters"+'\n')
input_file.write("    filename =  "+repr(filename)+'\n')
input_file.write("    ,fileout =  "+repr(fileout)+'\n')
input_file.write("    ,filename_noise =  "+repr(filename)+'\n')
input_file.write("    ,n_gauss =  "+repr(n_gauss)+'\n')
input_file.write("    ,lambda_amp =  "+repr(lambda_amp)+'d0'+'\n')
input_file.write("    ,lambda_mu =  "+repr(lambda_mu)+'d0'+'\n')
input_file.write("    ,lambda_sig =  "+repr(lambda_sig)+'d0'+'\n')
input_file.write("    ,lambda_var_amp =  "+repr(lambda_var_amp)+'d0'+'\n')
input_file.write("    ,lambda_var_mu =  "+repr(lambda_var_mu)+'d0'+'\n')
input_file.write("    ,lambda_var_sig =  "+repr(lambda_var_sig)+'d0'+'\n')
input_file.write("    ,amp_fact_init =  "+repr(amp_fact_init)+'d0'+'\n')
input_file.write("    ,sig_init =  "+repr(sig_init)+'d0'+'\n')
input_file.write("    ,init_option =  "+repr(init_option)+'\n')
input_file.write("    ,maxiter_init =  "+repr(maxiter_init)+'\n')
input_file.write("    ,maxiter =  "+repr(maxiter)+'\n')
input_file.write("    ,m =  "+repr(m)+'\n')
input_file.write("    ,noise =  "+noise+'\n')
input_file.write("    ,regul =  "+regul+'\n')
input_file.write("    ,descent =  "+descent+'\n')
input_file.write("    ,lstd =  "+repr(lstd)+'\n')
input_file.write("    ,ustd =  "+repr(ustd)+'\n')
input_file.write("    ,iprint =  "+repr(iprint)+'\n')
input_file.write("    ,iprint_init =  "+repr(iprint_init)+'\n')
input_file.write("    /"+'\n')
input_file.close()
# os.system("make run")
