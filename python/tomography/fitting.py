import sys, os
import numpy as np

from scipy.special import xlogy #, loggamma
from scipy.interpolate import PchipInterpolator

from tomography import *
from ..constants import *
from ..physics import *


###############################
#   Various fitting metrics   #
###############################

# Fitting metric for reconstructed Gamma bars
def fitting_metric_Gamma_bar(x, ebins, rbins, obs_gb, obs_gb_err):
    n_rbins = len(rbins)-1
    theo_gb = parametric_Gamma_bar(rbins, ebins, x[0], x[1:(1+n_rbins)], x[(1+n_rbins):])
    rdiff = (theo_gb[:,1:] - obs_gb[:,1:])/obs_gb_err[:,1:]
    return np.sum(rdiff*rdiff)

# Fitting metric for general polynomials
def fitting_metric_poly(x, rbins, matrix, p_moments, counts, exp_settings=iaxo):
    x = np.array(x)
    if np.any(x < 0):
        print("Problematic x =", x)
        raise ValueError("Some values of parameter point x are negative or NAN.")
    g = x[0]
    g2 = g*g 
    n_rbins = len(rbins)-1
    ch2 = 0 
    for p,o in zip(p_moments,counts):
        # Compute the values of the GBPs at rbins positions, and their max value
        ys = np.array([g2*p(x[j+1],x[n_rbins+j+1]) for j in range(n_rbins)]+[0])
        ys_max = np.max(ys)
        # Compute the polynomial coefficients for the rescaled ys (for numerical stability)
        try:
            coeffs = np.flip(PchipInterpolator(rbins, ys/ys_max).c, axis=0).reshape(-1)
        except ValueError:
            print("Problematic Gammabar detected: ", ys, flush=True)
            # Set spline to zero; TODO: could return np.inf or np.nan if opt algo accepts it
            coeffs = np.flip(PchipInterpolator(rbins, (len(rbins)+1)*[0]).c, axis=0).reshape(-1)
        t = ys_max*distance_factor*flux_to_events(matrix@coeffs, gagg=g, exp_settings=exp_settings)
        if np.any(t < 0):
            print("t =", t, " at x =", x)
            raise ValueError("Negative values for the rate coefficients occured.")
        ch2 += -2.0*np.sum(xlogy(o, t) - xlogy(o, o) - t + o)
    if (ch2 < 0):
        raise ValueError("Negative value for the fitting metric occured.")
    return ch2


#####################################
#  For the optimisation with emcee  #
#####################################

def log_prior(x, rbins, xb):
    n_rbinsp1 = len(rbins)
    truth_vals = [b[0] < x0 < b[1] for x0,b in zip(x,xb)]
    if sum(truth_vals) < len(truth_vals):
        return -np.inf
    xi = x[1:n_rbinsp1]/x[n_rbinsp1:]
    diff = (xi - xi_m)/xi_std
    return -0.5*np.sum(diff*diff)

def log_like(x, rbins, matrix, pm_interp, counts):
    return -0.5*fitting_metric_poly(x, rbins, matrix, pm_interp, counts)

def log_prob(x, rbins, counts, matrix, pm_interp, xb):
    lp = log_prior(x, rbins, xb)
    if not np.isfinite(lp):
        return -np.inf, -np.inf, np.nan
    ll = log_like(x, rbins, matrix, pm_interp, counts)
    return lp + ll, lp, ll