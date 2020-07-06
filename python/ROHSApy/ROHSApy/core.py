# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import ndimage
from scipy import optimize

plt.ion()
cm = plt.get_cmap('inferno')
cm.set_bad(color='black')
imkw = dict(origin='lower', interpolation='none', cmap=cm)

class ROHSA(object):
    def __init__(self, cube, hdr=None, filename=None):
        super(ROHSA, self).__init__()
        self.cube = cube
        self.hdr = hdr if hdr is not None else None
        self.filename = filename if filename is not None else "cube.fits"
        if self.hdr is not None : self.v = self.mean2vel(self.hdr["CRVAL3"]*1.e-3, self.hdr["CDELT3"]*1.e-3, self.hdr["CRPIX3"], np.arange(self.cube.shape[0]))


    def cube2dat(self, filename=None):
        if not filename : filename = self.filename[:-5] + ".dat" if self.filename is not None else "cube.dat"
        print("Generate " + filename + " file readable by fortran")
        
        with open(filename,'w+') as f:
            f.write('{:d}\t{:d}\t{:d}\n'.format(self.cube.shape[0], self.cube.shape[1], self.cube.shape[2]))
            for i in range(self.cube.shape[1]):
                for j in range(self.cube.shape[2]):
                    for k in range(self.cube.shape[0]):
                        line = '{:d}\t{:d}\t{:d}\t{:0.16f}\n'.format(k, i, j, self.cube[k,i,j])
                        f.write(line)


    def rms_map(self, rms_map, filename=None):
        if not filename :
            print("Generate rms_map.dat file")
        else: print("Generate " + filename + " file")

        filename = filename or "myrms_map.dat"
        
        with open(filename,'w+') as f:
            f.write('{:d}\t{:d}\n'.format(self.cube.shape[1], self.cube.shape[2]))
            for j in range(self.cube.shape[2]):
                for k in range(self.cube.shape[1]):
                    line = '{:d}\t{:d}\t{:0.16f}\n'.format(k, j, rms_map[k,j])
                    f.write(line)

        
    def gen_parameters(self, filename_parameters=None, filename=None, fileout="result.dat", timeout="timestep.dat", filename_noise="", n_gauss=3, lambda_amp=1000, 
                       lambda_mu=1000, lambda_sig=1000, lambda_var_amp=0, lambda_var_mu=0, lambda_var_sig=1000, lambda_lym_sig=0, amp_fact_init=0.66, sig_init = 4., 
                       lb_sig_init=1., ub_sig_init=100., lb_sig=0.001, ub_sig=100., init_option="mean", maxiter_init=15000, maxiter=800, m=10, 
                       noise=".false.", regul = ".true.", descent = ".true.", lstd = 1, ustd = 20, iprint = -1, iprint_init = -1, 
                       save_grid=".true.", lym=".false."):

        if not filename : 
            print("Need an input filename")
            sys.exit()

        if not filename_parameters :
            print("Generate parameters.txt file")
        else: print("Generate " + filename_parameters + " file")
        
        filename_parameters = filename_parameters or "parameters.txt"

        input_file = open(filename_parameters, 'w')
        input_file.write("&user_parameters"+'\n')
        input_file.write("    filename =  "+repr(filename)+'\n')
        input_file.write("    ,fileout =  "+repr(fileout)+'\n')
        input_file.write("    ,timeout =  "+repr(timeout)+'\n')
        input_file.write("    ,filename_noise =  "+repr(filename_noise)+'\n')
        input_file.write("    ,n_gauss =  "+repr(n_gauss)+'\n')
        input_file.write("    ,lambda_amp =  "+repr(lambda_amp)+'d0'+'\n')
        input_file.write("    ,lambda_mu =  "+repr(lambda_mu)+'d0'+'\n')
        input_file.write("    ,lambda_sig =  "+repr(lambda_sig)+'d0'+'\n')
        input_file.write("    ,lambda_var_amp =  "+repr(lambda_var_amp)+'d0'+'\n')
        input_file.write("    ,lambda_var_mu =  "+repr(lambda_var_mu)+'d0'+'\n')
        input_file.write("    ,lambda_var_sig =  "+repr(lambda_var_sig)+'d0'+'\n')
        input_file.write("    ,lambda_lym_sig =  "+repr(lambda_lym_sig)+'d0'+'\n')
        input_file.write("    ,amp_fact_init =  "+repr(amp_fact_init)+'d0'+'\n')
        input_file.write("    ,sig_init =  "+repr(sig_init)+'d0'+'\n')
        input_file.write("    ,lb_sig_init =  "+repr(lb_sig_init)+'d0'+'\n')
        input_file.write("    ,ub_sig_init =  "+repr(ub_sig_init)+'d0'+'\n')
        input_file.write("    ,lb_sig =  "+repr(lb_sig)+'d0'+'\n')
        input_file.write("    ,ub_sig =  "+repr(ub_sig)+'d0'+'\n')
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
        input_file.write("    ,lym =  "+lym+'\n')
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

    def physical_gaussian(self, gaussian):
        n_gauss = gaussian.shape[0]/3
        output = np.zeros(gaussian.shape)
        if self.hdr is not None :
            output[0::3] = gaussian[0::3]
            output[1::3] = self.mean2vel(self.hdr["CRVAL3"]*1.e-3, self.hdr["CDELT3"]*1.e-3, self.hdr["CRPIX3"], gaussian[1::3])
            output[2::3] = gaussian[2::3] * np.abs(self.hdr["CDELT3"])*1.e-3
            return output
        else:
            print("Missing header")
            return 0        

    def gauss(self, x, a, mu, sig):
        return a * np.exp(-((x - mu)**2)/(2. * sig**2))


    def gauss_2D(self, xs, a, mu, sig):
        return [a * np.exp(-((x - mu)**2)/(2. * sig**2)) for x in xs]

    
    def mean2vel(self, CRVAL, CDELT, CRPIX, mean):
        return [(CRVAL + CDELT * (mean[i] - CRPIX)) for i in range(len(mean))]            


    def plot_spect(self, gaussian, idy=0, idx=0):
        n_gauss = gaussian.shape[0]/3
        x = np.arange(self.cube.shape[0])

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        if self.hdr is not None :
            if not self.hdr["CRVAL3"] : print("Missing CRVAL3 keyword")
            v = self.mean2vel(self.hdr["CRVAL3"]*1.e-3, self.hdr["CDELT3"]*1.e-3, self.hdr["CRPIX3"], x)
            if v[0] > v[1] : v = v[::-1]
            ax.step(v, self.cube[:,idy,idx], color='cornflowerblue')
            tot = np.zeros(self.cube.shape[0])
            for i in np.arange(n_gauss):
                spectrum = self.gauss(x, gaussian[int(0+(3*i)),idy,idx], gaussian[int(1+(3*i)),idy,idx], gaussian[int(2+(3*i)),idy,idx])
                tot += spectrum
                ax.plot(v, spectrum, color="k")
            ax.plot(v, tot, color="r") 
            ax.set_ylabel(r'T [k]')
            ax.set_xlabel(r'v [km s$^{-1}$]')
        else:
            ax.step(x, self.cube[:,idy,idx], color='cornflowerblue')
            tot = np.zeros(self.cube.shape[0])
            for i in np.arange(n_gauss):
                spectrum = self.gauss(x, gaussian[int(0+(3*i)),idy,idx], gaussian[int(1+(3*i)),idy,idx], gaussian[int(2+(3*i)),idy,idx])
                tot += spectrum
                ax.plot(x, spectrum, color="k")
            ax.plot(x, tot, color="r") 
            ax.set_ylabel(r'T [k]')
            ax.set_xlabel(r'idx [pixel unit]')
             
        return 0 


    def return_result_cube(self, gaussian=None, ampfield=None, pixfield=None, sigfield=None):
        if gaussian is not None:
            result = np.zeros(self.cube.shape)
            n_gauss = gaussian.shape[0]/3
            for i in np.arange(n_gauss) :
                result += self.gauss_2D(np.arange(self.cube.shape[0]), gaussian[int(0+(3*i))], gaussian[int(1+(3*i))], gaussian[int(2+(3*i))])
            return result
        elif ampfield is not None and pixfield is not None and sigfield is not None:
            result = np.zeros(self.cube.shape)
            n_gauss = ampfield.shape[0]
            for i in np.arange(n_gauss) :
                result += self.gauss_2D(np.arange(self.cube.shape[0]), ampfield[i], pixfield[i], sigfield[i])
            return result
        else: print("error : 1 or 3 arguments needed")


    def write_fits(self, gaussian, fileout=None, hdr=False):
        if not fileout : fileout = self.filename[:-5] + "_ROHSA.fits" if self.filename is not None else "cube_ROHSA.fits"

        result = self.return_result_cube(gaussian)

        hdu = fits.PrimaryHDU(result)
        if self.hdr is not None : hdu.header = self.hdr
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(fileout, overwrite=True)        
        return 0
                
        
if __name__ == '__main__':    
    #Load data
    filename = "GHIGLS_DFN_Tb.fits"
    hdu = fits.open(filename)
    hdr = hdu[0].header
    cube = hdu[0].data[0][150:350,:4,:4]

    #Call ROHSApy
    core = ROHSA(cube, hdr=hdr, filename=filename)
    core.cube2dat()
    # core.rms_map_const()
    rms = np.std(cube[0:19,:,:], 0)
    core.rms_map(rms, filename="rms_map.dat")
    core.gen_parameters(filename="GHIGLS_DFN_Tb.dat", save_grid=".false.", iprint=1, 
                        lambda_var_sig=0, filename_noise="rms_map.dat", noise=".true.",
                        )
    core.run("parameters.txt", nohup=False)
    gaussian = core.read_gaussian("result.dat")
    # foo = core.write_fits(gaussian)
    # core.plot_spect(gaussian, 14, 14)

    def update_level(params, cube, std_map, n_gauss, lambda_amp, lambda_mu, lambda_sig, lambda_var_sig):
        bounds_inf = np.zeros(params.shape)
        bounds_sup = np.zeros(params.shape)
        
        for i in range(cube.shape[1]):
            for j in range(cube.shape[2]):
                for k in np.arange(params.shape[0]/3):
                    bounds_A = [0., np.max(cube[:,i,j])]
                    bounds_mu = [0., len(cube[:,i,j])]
                    bounds_sigma = [1., 100.]
                    
                    bounds_inf[int(0+(3*k)),i,j] = bounds_A[0]
                    bounds_inf[int(1+(3*k)),i,j] = bounds_mu[0]
                    bounds_inf[int(2+(3*k)),i,j] = bounds_sigma[0]
                    
                    bounds_sup[int(0+(3*k)),i,j] = bounds_A[1]
                    bounds_sup[int(1+(3*k)),i,j] = bounds_mu[1]
                    bounds_sup[int(2+(3*k)),i,j] = bounds_sigma[1]
                
        bounds = [(bounds_inf.ravel()[i], bounds_sup.ravel()[i]) for i in np.arange(len(bounds_sup.ravel()))]
        result = optimize.minimize(cost, params.ravel(), method="L-BFGS-B", jac=False, args=(cube, std_map, n_gauss, lambda_amp, 
                                                                                             lambda_mu, lambda_sig, lambda_var_sig), 
                                   bounds=bounds, options={'maxiter': 2})
        return result


    def gauss(x, A, mu, sigma):
        return A * np.exp(-((x - mu)**2)/(2. * sigma**2))


    def cost(pars, data, std_map, n_gauss, lambda_amp, lambda_mu, lambda_sig, lambda_var_sig):
        params = np.reshape(pars, (3*n_gauss, data.shape[1], data.shape[2]))
        
        x = np.arange(data.shape[0])
        
        model = np.zeros(data.shape)
        
        for i in np.arange(data.shape[1]):
            for j in np.arange(data.shape[2]):
                for k in np.arange(n_gauss):
                    line = gauss(x, params[0+(k*3),i,j], params[1+(k*3),i,j], params[2+(k*3),i,j])
                    model[:,i,j] += line

        J = np.sum((((model - data)/std_map).ravel())**2)
        
        kernel = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]]) / 4.
        
        R = 0.
        for k in np.arange(n_gauss):
            conv_amp = ndimage.convolve(params[0+(k*3),:,:], kernel, mode='reflect')
            conv_mu = ndimage.convolve(params[1+(k*3),:,:], kernel, mode='reflect')
            conv_sig = ndimage.convolve(params[2+(k*3),:,:], kernel, mode='reflect')

            R += lambda_amp * np.sum(conv_amp**2)
            R += lambda_mu * np.sum(conv_mu**2)
            R += lambda_sig * np.sum(conv_sig**2)
            # R += lambda_var_sig * np.sum((params[2+(k*3),:,:] - np.mean(params[2+(k*3),:,:]))**2)
            
        return 0.5*J + 0.5*R

    result = update_level(gaussian, cube, rms, 3, 1000,1000, 1000, 0)
    hess_inv = result.hess_inv.todense()
    err = np.sqrt(np.diag(hess_inv)) / cost(gaussian, cube, std_map, 3, 1000,1000, 1000., 0)  
