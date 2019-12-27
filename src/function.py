import numpy as np
from astropy import constants as const
from astropy import units as u

cexp = (const.h.cgs*const.c.cgs/const.k_B.cgs).value

def planck_l(l, td):
    return 1./(np.exp(cexp/l/td) - 1.)

def MBB_l(l,tau,beta,td,l0):
    return tau * (l0/l)**beta * 1./(np.exp(cexp/l/td) - 1.) #planck_l(l,td)

def d_MBB_l_d_tau(l,tau,beta,td,l0):
    return (l0/l)**beta * planck_l(l,td)

def d_MBB_l_d_beta(l,tau,beta,td,l0):
    return tau * np.log(l0/l) * (l0/l)**beta * planck_l(l,td)

def d_MBB_l_d_td(l,tau,beta,td,l0):  
    a = cexp / l
    a_td = a / td
    return tau * (l0/l)**beta * a * np.exp(a_td) / (td**2. * (np.exp(a_td) - 1.)**2.)

def butterworth(k, H0, k0, n):
    return H0 / np.sqrt((1. + (k/k0)**(2*n)))

def poly(coef, xx, yy, degree):
    model = 0.
    k = 0
    for i in np.arange(degree):
        for j in np.arange(degree+1-i):
            model += coef[k] * xx**i * yy**j
            k += 1
    return model

def d_poly_d_xx(coef, xx, yy, degree):
    model = 0.
    k = 0
    for i in np.arange(degree):
        for j in np.arange(degree+1-i):
            model += coef[k] * i * xx**(i-1) * yy**j
            k += 1
    return model

def d_poly_d_yy(coef, xx, yy, degree):
    model = 0.
    k = 0
    for i in np.arange(degree):
        for j in np.arange(degree+1-i):
            model += coef[k] * j * xx**i * yy**(j-1)
            k += 1
    return model

def colorcorr(color,l,beta,td,degree):
    return np.array([poly(color[i,:],beta,td,degree) for i in np.arange(len(color))])

def d_colorcorr_d_beta(color,l,beta,td,degree):
    return np.array([d_poly_d_xx(color[i,:],beta,td,degree) for i in np.arange(len(color))])

def d_colorcorr_d_td(color,l,beta,td,degree):
    return np.array([d_poly_d_yy(color[i,:],beta,td,degree) for i in np.arange(len(color))])
    
def MBB_l_cc(l,tau,beta,td,l0,color,degree):    
    cc = colorcorr(color,l,beta,td,degree) 
    mbb = MBB_l(l,tau,beta,td,l0)
    return  np.array([cc[i] * mbb[i] for i in np.arange(len(color))])

def d_MBB_l_cc_d_tau(l,tau,beta,td,l0,color,degree):    
    cc = colorcorr(color,l,beta,td,degree) 
    mbb = MBB_l(l,tau,beta,td,l0)
    d_mbb_d_tau = d_MBB_l_d_tau(l,tau,beta,td,l0)    
    return  np.array([cc[i]*d_mbb_d_tau[i] for i in np.arange(len(color))])

def d_MBB_l_cc_d_beta(l,tau,beta,td,l0,color,degree):    
    cc = colorcorr(color,l,beta,td,degree) 
    d_color_d_beta = d_colorcorr_d_beta(color,l,beta,td,degree) 
    mbb = MBB_l(l,tau,beta,td,l0)
    d_mbb_d_beta = d_MBB_l_d_beta(l,tau,beta,td,l0)    
    return  np.array([d_color_d_beta[i] * mbb[i] + cc[i] * d_mbb_d_beta[i] for i in np.arange(len(color))])

def d_MBB_l_cc_d_td(l,tau,beta,td,l0,color,degree):    
    cc = colorcorr(color,l,beta,td,degree) 
    d_cc_d_td = d_colorcorr_d_td(color,l,beta,td,degree) 
    mbb = MBB_l(l,tau,beta,td,l0)
    d_mbb_l_d_td = d_MBB_l_d_td(l,tau,beta,td,l0)    
    return  np.array([d_cc_d_td[i] * mbb[i] + cc[i] * d_mbb_l_d_td[i] for i in np.arange(len(color))])

