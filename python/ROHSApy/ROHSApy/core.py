# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

plt.ion()
cm = plt.get_cmap('inferno')
cm.set_bad(color='black')
imkw = dict(origin='lower', interpolation='none', cmap=cm)

class ROHSA(object):
    def __init__(self, cube, filename=None):
        super(ROHSA, self).__init__()
        self.cube = cube


    def cube2dat(self, filename=None):
        if not filename :
            print("Generate mycube.dat file")
        else: print("Generate " + filename + " file")

        filename = filename or "mycube.dat"
        
        with open(filename,'w+') as f:
            f.write('{:d}\t{:d}\t{:d}\n'.format(self.cube.shape[0], self.cube.shape[1], self.cube.shape[2]))
            for i in range(self.cube.shape[1]):
                for j in range(self.cube.shape[2]):
                    for k in range(self.cube.shape[0]):
                        line = '{:d}\t{:d}\t{:d}\t{:0.16f}\n'.format(k, i, j, self.cube[k,i,j])
                        f.write(line)


    def rms_map_const(self, filename=None, const=1.):
        if not filename :
            print("Generate rms_map.dat file")
        else: print("Generate " + filename + " file")

        filename = filename or "myrms_map.dat"
        
        with open(filename,'w+') as f:
            f.write('{:d}\t{:d}\n'.format(self.cube.shape[1], self.cube.shape[2]))
            for j in range(self.cube.shape[2]):
                for k in range(self.cube.shape[1]):
                    line = '{:d}\t{:d}\t{:0.16f}\n'.format(j, k, const)
                    f.write(line)

        
    def gen_parameters(self, filename_parameters=None, filename=None, fileout="result.dat", filename_noise="", n_gauss=3, lambda_amp=1000, 
                       lambda_mu=1000, lambda_sig=1000, lambda_var_amp=0, lambda_var_mu=0, lambda_var_sig=1000, amp_fact_init=0.66, sig_init = 4., 
                       init_option="mean", maxiter_init=15000, maxiter=800, m=10, noise=".false.", regul = ".true.", descent = ".true.",
                       lstd = 1, ustd = 20, iprint = -1, iprint_init = -1, save_grid=".true."):

        if not filename : 
            print("Need an input filename")
            sys.exit()

        if not filename_parameters :
            print("Generate parameters.dat file")
        else: print("Generate " + filename_parameters + " file")
        
        filename_parameters = filename_parameters or "parameters.txt"

        input_file = open(filename_parameters, 'w')
        input_file.write("&user_parameters"+'\n')
        input_file.write("    filename =  "+repr(filename)+'\n')
        input_file.write("    ,fileout =  "+repr(fileout)+'\n')
        input_file.write("    ,filename_noise =  "+repr(filename_noise)+'\n')
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
        input_file.write("    ,save_grid =  "+save_grid+'\n')
        input_file.write("    /"+'\n')
        input_file.close()


    def run(self, filename=None, nohup=False):
        if not filename: 
            print("Need the input filename parameters to run ROHSA")
            sys.exit()
        if nohup == False:
            os.system("ROHSA " + filename)
        else:
            os.system("nohup ROHSA " + filename + "&")


    def read_gaussian(self, filename=None):
        if not filename: 
            print("Need the ouput filename of ROHSA")
            sys.exit()

        data = np.genfromtxt(filename)
        print("Opening data file")

        amp = data[:, 2]
        mean = data[:, 3] - 1
        sigma = data[:, 4]
        
        dim_y = self.cube.shape[1]
        dim_x = self.cube.shape[2]

        n_gauss = int(len(amp) / (dim_y*dim_x))
        params = np.zeros((3*n_gauss, dim_y, dim_x))
        
        i__ = 0
        for i in range(dim_y):
            for j in range(dim_x):
                for k in range(n_gauss):
                    params[0+(3*k),i,j] = amp[i__]
                    params[1+(3*k),i,j] = mean[i__]
                    params[2+(3*k),i,j] = sigma[i__]
                    i__ += 1
        
        return params


    def gauss(self, x, a, mu, sig):
        return a * np.exp(-((x - mu)**2)/(2. * sig**2))


    def plot_spect(self, gaussian, idy=0, idx=0):
        plt.figure()
        x = np.arange(self.cube.shape[0])
        plt.plot(x, self.cube[:,idy,idx])
        tot = np.zeros(self.cube.shape[0])
        for i in np.arange(gaussian.shape[0]/3):
            spectrum = self.gauss(x, gaussian[int(0+(3*i)),idy,idx], gaussian[int(1+(3*i)),idy,idx], gaussian[int(2+(3*i)),idy,idx])
            tot += spectrum
            plt.plot(x, spectrum, color="k") 
        plt.plot(x, tot, color="r") 
             
        return 0 
                
        
if __name__ == '__main__':    
    #Load data
    filename = "GHIGLS_DFN_Tb.fits"
    hdu = fits.open(filename)
    cube = hdu[0].data[0][150:350,:32,:32]

    #Call ROHSApy
    core = ROHSA(cube)
    core.cube2dat()
    core.rms_map_const()
    core.gen_parameters(filename="mycube.dat", save_grid=".false.")
    core.run("parameters.txt", nohup=False)
    gaussian = core.read_gaussian("result.dat")
    core.plot_spect(gaussian, 14, 14)
