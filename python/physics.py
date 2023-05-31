import numpy as np

from scipy.special import exprel, xlog1py
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad

from constants import *

################################################
#  Parametric axion production and conversion  #
################################################

# Primakoff rate
def parametric_primakoff_rate(om, gagg, ks, temp):
    if (ks > 0) & (temp > 0) & (om > 0):
        x = om/temp
        # x = 0.5*om/temp
        y = 1e-16*gagg*ks*temp
        z = 0.5*ks/om
        z2 = z*z
        w = exprel(x)
        # w = np.exp(-x)/np.real(np.sinc(x*1.0j/np.pi))
        # w = np.exp(-x)
        # w = w/(1.0 - w)
        # return y*y*(xlog1py(1.0+z2, 1.0/z2) - 1.0)*w/(16.0*np.pi*om)
        # return y*y*(xlog1py(1.0+z2, 1.0/z2) - 1.0)*w/(16.0*np.pi*temp)
        bracket = xlog1py(1.0+z2, 1.0/z2) - 1.0
        denom = 16.0*np.pi*om*w
        return (y*bracket)*(y/denom)
    else:
        return 0

# Conversion probability
def p_conv_m0(gagg, bfield, length):
    x = 1e-16*gagg*bfield*length
    return 0.25*x*(x*c_p_conv)

# Calculate the expected number of events in IAXO from the flux
def flux_to_events(grid, gagg=0.5, exp_settings=iaxo):
    p_conv = p_conv_m0(gagg, exp_settings['B'], exp_settings['L'])
    flux_to_events = exp_settings['eff']*exp_settings['A']*exp_settings['t']*p_conv
    return grid*flux_to_events


#############################################
#  Pre-computed Primakoff moments in ebins  #
#############################################

def primakoff_moments_integral(ebins):
    n_ebins = len(ebins)-1
    pm_full = [lambda k, t, j=j: quad(lambda om: om*om*parametric_primakoff_rate(om, 1, k, t), ebins[j], ebins[j+1])[0]/(2*np.pi*np.pi) for j in range(n_ebins)]
    return pm_full

def primakoff_moments(ebins, nk=2001, nt=1001):
    n_ebins = len(ebins)-1
    kvals = np.linspace(0, 20, nk)
    tvals = np.linspace(0, 4, nt)
    pintegral = lambda k, t, e1, e2: quad(lambda om: om*om*parametric_primakoff_rate(om, 1, k, t), e1, e2)[0]/(2*np.pi*np.pi)
    primakoff_moments = np.array([[[pintegral(k, t, ebins[i], ebins[i+1]) for t in tvals] for k in kvals] for i in range(n_ebins)])
    pm_interp = [RegularGridInterpolator((kvals, tvals), p, method='linear') for p in primakoff_moments]
    pm_interp = [lambda k, t, p=p: p((k,t)) for p in pm_interp]
    return pm_interp