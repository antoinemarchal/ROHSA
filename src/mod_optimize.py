import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy import ndimage
from scipy import optimize


def init_reshape(cube, nside):
    """! 
    Reshape cube around the center of the initial map 
    
    @param cube: 3D array : hyperspectral data (PPV)
    @param nside: Int : Nside of map  / \f$ 2^{Nside} \f$

    @return reshape cube
    """
    half_nside = 2**nside / 2
    
    dim1 = cube.shape[1]/2 - half_nside
    dim2 = cube.shape[1]/2 + half_nside
    dim3 = cube.shape[2]/2 - half_nside
    dim4 = cube.shape[2]/2 + half_nside

    return cube[:, dim1:dim2, dim3:dim4]


def go_up_level(cube_params):
    dim = cube_params.shape
    cube_params_up = np.zeros((dim[0],dim[1]*2,dim[2]*2))
    for i in range(cube_params.shape[1]):
        for j in range(cube_params.shape[2]):
            for k in range(2):
                for l in range(2):
                    cube_params_up[:,k+((i)*2),l+((j)*2)] = cube_params[:,i,j]
    return cube_params_up


def mean_cube(cube,nside):
    spectrum = 0.
    n = cube.shape[1] / nside
    cube_mean = np.zeros((cube.shape[0], nside, nside))
    for i in range(cube_mean.shape[1]):
        for j in range(cube_mean.shape[2]):
            for k in range(n):
                for l in range(n):
                    spectrum = spectrum + cube[:,k+((i)*n),l+((j)*n)]
            spectrum = spectrum / (n**2)
            cube_mean[:,i,j] = spectrum
            spectrum = 0.

    return cube_mean


def MBB_nu(nu,NHI,sigma,beta,T,nu0):
    return sigma * (nu/nu0)**beta * NHI * planck_nu(nu,T)


def MBB_l(l,NHI,sigma,beta,T,l0):
    return sigma * (l/l0)**beta * NHI * planck_l(l,T)


def MBB_nu_adim(nu,NHI,sigma,beta,T,nu0):
    return sigma * (nu/nu0)**beta * NHI * planck_nu(nu,T) / (2. * const.h * nu**3. / const.c**2)


def MBB_l_adim(l,NHI,sigma,beta,T,l0):
    return sigma * (l0/l)**beta * NHI * planck_l_adim(l,T)


def planck_nu(nu, T):
    return 2. * const.h * nu**3. / const.c**2 * 1./(np.exp(const.h*nu/const.k_B/T) - 1.)


def planck_l(l, T):
    return 2. * const.h * const.c**2 / l**5. * 1./(np.exp(const.h*const.c/l/const.k_B/T) - 1.)


def planck_l_adim(l, T):
    return 1./(np.exp(const.h*const.c/l/const.k_B/T) - 1.)


def f_g_mean(params, n_mbb, wave, data, l0, NHI):
    model = np.zeros(data.shape)
    F = np.zeros(data.shape)
    dF_over_dB = np.zeros((len(params), len(data)))
    product = np.zeros((len(params), len(data)))
    deriv = np.zeros((len(params)))

    for k in np.arange(n_mbb):
        line = MBB_l_adim(wave,NHI[k],params[0+(k*3)]*u.cm**2, params[1+(k*3)], params[2+(k*3)]*u.K,l0)
        model += line.value
    
        dF_over_dB[0+(k*3)] = ((l0/wave)**params[1+(k*3)] * NHI[k] * planck_l_adim(wave,params[2+(k*3)]*u.K)).to(u.cm**-2)
        dF_over_dB[1+(k*3)] = params[0+(k*3)]*u.cm**2 * np.log(l0/wave) * (l0/wave)**params[1+(k*3)] * NHI[k] * planck_l_adim(wave,params[2+(k*3)]*u.K)
        dF_over_dB[2+(k*3)] = (params[0+(k*3)]*u.cm**2 * (l0/wave)**params[1+(k*3)] * NHI[k] * (const.h*const.c/wave/const.k_B) * np.exp(const.h*const.c/wave/const.k_B/(params[2+(k*3)]*u.K)) * 1./ ((params[2+(k*3)]*u.K)**2. * (np.exp(const.h*const.c/wave/const.k_B/(params[2+(k*3)]*u.K)) - 1.)**2.)).to(u.K**-1)

    J = np.sum(((model - data).ravel())**2)
    F = model - data

    product = dF_over_dB * F
    deriv = np.sum(product, axis=1)    
    
    return 0.5*J, deriv


def f_g(pars, n_mbb, wave, data, l0, NHI, lambda_sig, lambda_beta, lambda_T, lambda_var_sig, lambda_var_beta, lambda_var_T, kernel): 
    dim_params = n_mbb*3*data.shape[1]*data.shape[2]
    params = np.reshape(pars[:dim_params], (3*n_mbb, data.shape[1], data.shape[2]))
    b = np.zeros(n_mbb)
    c = np.zeros(n_mbb)
    d = np.zeros(n_mbb)
    
    for i in np.arange(len(b)): 
        b[i] = pars[dim_params+(0+(3*i))]
        c[i] = pars[dim_params+(1+(3*i))]
        d[i] = pars[dim_params+(2+(3*i))]
        
    # print np.mean(params[0::3],(1,2))
    # print np.mean(params[1::3],(1,2))
    # print np.mean(params[2::3],(1,2))
    
    model = np.zeros(data.shape)
    F = np.zeros(data.shape)
    dF_over_dB = np.zeros((params.shape[0], data.shape[0], data.shape[1], data.shape[2]))
    conv = np.zeros((params.shape[0], data.shape[1], data.shape[2]))
    dR_over_dB = np.zeros((params.shape[0], data.shape[1], data.shape[2]))
    grad_cube = np.zeros((params.shape[0], data.shape[1], data.shape[2]))
    R = 0.
    R_var = 0.
    
    for k in np.arange(n_mbb):
        line = MBB_l_adim(wave,NHI[k],params[0+(k*3)]*u.cm**2, params[1+(k*3)], params[2+(k*3)]*u.K,l0)
        model += line.value
    
        dF_over_dB[0+(k*3)] = ((l0/wave)**params[1+(k*3)] * NHI[k] * planck_l_adim(wave,params[2+(k*3)]*u.K)).to(u.cm**-2)
        dF_over_dB[1+(k*3)] = params[0+(k*3)]*u.cm**2 * np.log(l0/wave) * (l0/wave)**params[1+(k*3)] * NHI[k] * planck_l_adim(wave,params[2+(k*3)]*u.K)
        dF_over_dB[2+(k*3)] = (params[0+(k*3)]*u.cm**2 * (l0/wave)**params[1+(k*3)] * NHI[k] * (const.h*const.c/wave/const.k_B) * np.exp(const.h*const.c/wave/const.k_B/(params[2+(k*3)]*u.K)) * 1./ ((params[2+(k*3)]*u.K)**2. * (np.exp(const.h*const.c/wave/const.k_B/(params[2+(k*3)]*u.K)) - 1.)**2.)).to(u.K**-1)

        conv[0+(k*3)] = lambda_sig * ndimage.convolve(params[0+(k*3)], kernel, mode='reflect')
        conv[1+(k*3)] = lambda_beta * ndimage.convolve(params[1+(k*3)], kernel, mode='reflect')
        conv[2+(k*3)] = lambda_T * ndimage.convolve(params[2+(k*3)], kernel, mode='reflect')

        dR_over_dB[0+(k*3)] = lambda_sig * ndimage.convolve(conv[0+(k*3)], kernel, mode='reflect')
        dR_over_dB[1+(k*3)] = lambda_beta * ndimage.convolve(conv[1+(k*3)], kernel, mode='reflect')
        dR_over_dB[2+(k*3)] = lambda_T * ndimage.convolve(conv[2+(k*3)], kernel, mode='reflect')

    J = np.sum(((model - data).ravel())**2)
    R = np.sum(conv**2)    

    for k in np.arange(len(b)): 
        R_var += lambda_var_sig * np.sum((params[0+(k*3)] - b[k])**2.)
        R_var += lambda_var_beta * np.sum((params[1+(k*3)] - c[k])**2.)
        R_var += lambda_var_T * np.sum((params[2+(k*3)] - d[k])**2.)

    F = model - data

    grad = np.zeros(pars.shape)
    grad_F_times_F = np.sum(dF_over_dB * F,axis=1)
    grad_R = dR_over_dB
    grad_R_var_sig = (params[0::3].transpose() - b).transpose()
    grad_R_var_beta = (params[1::3].transpose() - c).transpose()
    grad_R_var_T = (params[2::3].transpose() - d).transpose()

    grad_cube = grad_F_times_F + grad_R
    grad_cube[0::3] += lambda_var_sig * grad_R_var_sig
    grad_cube[1::3] += lambda_var_beta * grad_R_var_beta
    grad_cube[2::3] += lambda_var_T * grad_R_var_T

    grad[:dim_params] = grad_cube.ravel()
    for k in np.arange(len(b)): 
        grad[dim_params+(0+(k*3))] = - lambda_var_sig * np.sum((params[0+(k*3)] - b[k]))
        grad[dim_params+(1+(k*3))] = - lambda_var_beta * np.sum((params[1+(k*3)] - c[k]))
        grad[dim_params+(2+(k*3))] = - lambda_var_T * np.sum((params[2+(k*3)] - d[k]))
        
    return 0.5*(J + R + R_var), grad


def update_level(nside, cube, params, wave, NHI, l0, lambda_sig, lambda_beta, lambda_T, bounds, maxiter, lb_sig, ub_sig, lb_beta, ub_beta, lb_T, ub_T):
    bounds = init_bounds(cube, params, lb_sig, ub_sig, lb_beta, ub_beta, lb_T, ub_T)

    result = optimize.fmin_l_bfgs_b(f_g, params.ravel(), args=(n_mbb, wave, cube, l0, NHI, lambda_sig, lambda_beta, lambda_T), 
                                    bounds=bounds, approx_grad=False, disp=1, maxiter=maxiter)

    return np.reshape(result[0], (params.shape))

    
def init_bounds(cube, params, lb_sig, ub_sig, lb_beta, ub_beta, lb_T, ub_T):
    bounds_inf = np.zeros(params.shape)
    bounds_sup = np.zeros(params.shape)

    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            for k in np.arange(params.shape[0]/3):
                bounds_A = [lb_sig, ub_sig]
                bounds_mu = [lb_beta, ub_beta]
                bounds_sigma = [lb_T, ub_T]
                
                bounds_inf[0+(3*k),i,j] = bounds_A[0]
                bounds_inf[1+(3*k),i,j] = bounds_mu[0]
                bounds_inf[2+(3*k),i,j] = bounds_sigma[0]
                
                bounds_sup[0+(3*k),i,j] = bounds_A[1]
                bounds_sup[1+(3*k),i,j] = bounds_mu[1]
                bounds_sup[2+(3*k),i,j] = bounds_sigma[1]
            
    return [(bounds_inf.ravel()[i], bounds_sup.ravel()[i]) for i in np.arange(len(bounds_sup.ravel()))]


def init_bounds_mean(n_mbb, lb_sig, ub_sig, lb_beta, ub_beta, lb_T, ub_T):
    bounds_inf = np.zeros(3*n_mbb)
    bounds_sup = np.zeros(3*n_mbb)
    
    bounds_inf[0::3] = lb_sig
    bounds_inf[1::3] = lb_beta
    bounds_inf[2::3] = lb_T
    
    bounds_sup[0::3] = ub_sig
    bounds_sup[1::3] = ub_beta
    bounds_sup[2::3] = ub_T
    
    return [(bounds_inf[i], bounds_sup[i]) for i in np.arange(len(bounds_sup))]


def reshape_cube_up(cube, nside):
    subcube = np.zeros((cube.shape[0], 2**nside, 2**nside))
    cube_w, cube_h = cube.shape[1], cube.shape[2]
    subcube_w, subcube_h = 2**nside, 2**nside
    offset = ((subcube_w - cube_w) / 2, (subcube_h - cube_h) / 2)
    
    subcube[:, offset[0]:offset[0]+cube_w, offset[1]:offset[1]+cube_h] = cube

    return subcube


# def reshape_cube_down(subcube, shape, nside): Je sais pas si ca marche
#     cube = np.zeros(shape)
#     cube_w, cube_h = cube.shape[1], cube.shape[2]
#     subcube_w, subcube_h = 2**nside, 2**nside
#     offset = ((subcube_w - cube_w) / 2, (subcube_h - cube_h) / 2)
    
#     cube = subcube[:, offset[0]:offset[0]+cube_w, offset[1]:offset[1]+cube_h]

#     return cube



